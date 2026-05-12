// Single-core ONNX Runtime C++ benchmark for the shift+smear-reweight
// polyhead ONNX produced by ``shift_smear_reweight_export.py``.
//
// Same benchmark structure as bench_aoti_shift_smear_reweight.cpp but
// uses ORT's CPU EP. Times trunk(y, c) -> coefs and the C++
// polynomial evaluation at N (u, σ) values per event, both
// individually and combined.
//
// Inputs:
//     <model.onnx> <indices.json>  [--batch B]  [--n-var N]
//                                  [--duration S]  [--cpu-pin C]
//
// Build:
//   g++ -O2 -std=c++17 \
//       -I /usr/include/onnxruntime \
//       -I scripts/corrections/muon_calibration \
//       bench_ort_shift_smear_reweight.cpp \
//       -lonnxruntime -lpthread \
//       -o bench_ort_shift_smear_reweight

#define _GNU_SOURCE
#include <sched.h>

#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include <onnxruntime_cxx_api.h>

#include "evaluate_joint.hpp"

static void pin_to_cpu(int cpu) {
    if (cpu < 0) return;
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(cpu, &mask);
    if (sched_setaffinity(0, sizeof(mask), &mask) != 0) {
        std::fprintf(stderr, "warning: sched_setaffinity(%d) failed\n", cpu);
    }
}

int main(int argc, char** argv) try {
    if (argc < 3) {
        std::fprintf(stderr,
            "usage: %s <model.onnx> <indices.json> "
            "[--batch B] [--n-var N] [--duration S] [--cpu-pin C]\n",
            argv[0]);
        return 2;
    }
    const char* onnx_path = argv[1];
    const char* json_path = argv[2];
    int64_t B     = 1;
    int64_t N_var = 32;
    double duration = 2.0;
    int cpu_pin = 0;
    for (int i = 3; i < argc; ++i) {
        if (std::strcmp(argv[i], "--batch") == 0 && i + 1 < argc) {
            B = std::atoll(argv[++i]);
        } else if (std::strcmp(argv[i], "--n-var") == 0 && i + 1 < argc) {
            N_var = std::atoll(argv[++i]);
        } else if (std::strcmp(argv[i], "--duration") == 0 && i + 1 < argc) {
            duration = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "--cpu-pin") == 0 && i + 1 < argc) {
            cpu_pin = std::atoi(argv[++i]);
        }
    }

    setenv("OMP_NUM_THREADS",      "1", 1);
    setenv("MKL_NUM_THREADS",      "1", 1);
    setenv("OPENBLAS_NUM_THREADS", "1", 1);
    pin_to_cpu(cpu_pin);

    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "bench_ssr");
    Ort::SessionOptions opts;
    opts.SetIntraOpNumThreads(1);
    opts.SetInterOpNumThreads(1);
    opts.SetExecutionMode(ORT_SEQUENTIAL);
    opts.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    opts.EnableCpuMemArena();

    Ort::Session session(env, onnx_path, opts);
    Ort::AllocatorWithDefaultOptions alloc;

    const size_t n_inputs  = session.GetInputCount();
    const size_t n_outputs = session.GetOutputCount();
    if (n_inputs != 2 || n_outputs < 1) {
        std::fprintf(stderr, "expected 2 inputs / >=1 output; got %zu/%zu\n",
                     n_inputs, n_outputs);
        return 3;
    }
    std::vector<std::string> in_names, out_names;
    for (size_t i = 0; i < n_inputs; ++i) {
        auto n = session.GetInputNameAllocated(i, alloc);
        in_names.emplace_back(n.get());
    }
    for (size_t i = 0; i < n_outputs; ++i) {
        auto n = session.GetOutputNameAllocated(i, alloc);
        out_names.emplace_back(n.get());
    }

    auto spec = ssr::load_basis_spec(json_path);

    auto infer_dim1 = [&](size_t idx, int64_t fallback) -> int64_t {
        Ort::TypeInfo ti = session.GetInputTypeInfo(idx);
        auto info = ti.GetTensorTypeAndShapeInfo();
        const size_t rank = info.GetDimensionsCount();
        if (rank < 2) return fallback;
        std::vector<int64_t> shape(rank);
        info.GetDimensions(shape.data(), rank);
        return shape[1] > 0 ? shape[1] : fallback;
    };
    const int64_t n_y = infer_dim1(0, spec.n_features);
    const int64_t n_c = infer_dim1(1, 5);
    std::fprintf(stderr,
        "model: %s\ninputs: %s[N,%ld], %s[N,%ld]  n_basis=%d\n"
        "B=%ld  N_var=%ld\n",
        onnx_path,
        in_names[0].c_str(), (long)n_y,
        in_names[1].c_str(), (long)n_c,
        spec.n_basis, (long)B, (long)N_var);

    auto mem_info = Ort::MemoryInfo::CreateCpu(
        OrtArenaAllocator, OrtMemTypeDefault);

    std::mt19937 rng(0);
    std::normal_distribution<float> gauss(0.0f, 1.0f);

    std::vector<float> y_buf(B * n_y), c_buf(B * n_c);
    for (auto& v : y_buf) v = gauss(rng);
    for (auto& v : c_buf) v = gauss(rng);

    std::array<int64_t, 2> y_shape{B, n_y};
    std::array<int64_t, 2> c_shape{B, n_c};

    Ort::Value y_tensor = Ort::Value::CreateTensor<float>(
        mem_info, y_buf.data(), y_buf.size(),
        y_shape.data(), y_shape.size());
    Ort::Value c_tensor = Ort::Value::CreateTensor<float>(
        mem_info, c_buf.data(), c_buf.size(),
        c_shape.data(), c_shape.size());

    // Bind output to a host-owned buffer so we can read coefs from
    // the same buffer iteration after iteration without re-allocating
    // (matches the AOTI bench shape).
    std::vector<float> coefs_buf((std::size_t)B * spec.n_basis);
    std::array<int64_t, 2> coefs_shape{B, spec.n_basis};
    Ort::Value coefs_tensor = Ort::Value::CreateTensor<float>(
        mem_info, coefs_buf.data(), coefs_buf.size(),
        coefs_shape.data(), coefs_shape.size());

    Ort::IoBinding binding(session);
    binding.BindInput(in_names[0].c_str(), y_tensor);
    binding.BindInput(in_names[1].c_str(), c_tensor);
    binding.BindOutput(out_names[0].c_str(), coefs_tensor);

    // Polynomial-side scratch.
    std::uniform_real_distribution<float> uni(-0.7f, 0.7f);
    std::vector<float> u_vars((std::size_t)B * N_var * n_y);
    std::vector<float> s_vars((std::size_t)B * N_var * n_y);
    std::vector<float> d_out ((std::size_t)B * N_var);
    for (auto& v : u_vars) v = uni(rng);
    for (auto& v : s_vars) v = uni(rng);

    auto run_loop = [&](double dur, auto&& body) {
        const int64_t min_iters = 5;
        int64_t n_iters = 0;
        auto t0 = std::chrono::steady_clock::now();
        auto deadline = t0 + std::chrono::duration_cast<
            std::chrono::steady_clock::duration>(
            std::chrono::duration<double>(dur));
        while (true) {
            body();
            ++n_iters;
            if (n_iters >= min_iters &&
                std::chrono::steady_clock::now() >= deadline) break;
        }
        auto t1 = std::chrono::steady_clock::now();
        return std::pair<int64_t, double>(
            n_iters,
            std::chrono::duration<double>(t1 - t0).count());
    };

    // Warmup.
    for (int i = 0; i < 50; ++i) {
        session.Run(Ort::RunOptions{nullptr}, binding);
        ssr::evaluate_joint_batched(
            spec, coefs_buf.data(),
            u_vars.data(), s_vars.data(),
            (int)B, (int)N_var, d_out.data());
    }

    // Trunk only.
    {
        auto [n, t] = run_loop(duration, [&]() {
            session.Run(Ort::RunOptions{nullptr}, binding);
        });
        double ms_b = 1000.0 * t / (double)n;
        double us_e = 1e6 * t / (double)n / (double)B;
        std::printf(
            "[trunk-only ] B=%-4ld iters=%-7ld %8.4f ms/call  "
            "%8.3f us/event\n",
            (long)B, (long)n, ms_b, us_e);
    }
    // Refresh coefs once so poly-only works on a real coefs vector.
    session.Run(Ort::RunOptions{nullptr}, binding);

    // Poly only.
    {
        auto [n, t] = run_loop(duration, [&]() {
            ssr::evaluate_joint_batched(
                spec, coefs_buf.data(),
                u_vars.data(), s_vars.data(),
                (int)B, (int)N_var, d_out.data());
        });
        double ms_b   = 1000.0 * t / (double)n;
        double us_e   = 1e6   * t / (double)n / (double)B;
        double us_var = us_e / (double)N_var;
        std::printf(
            "[poly-only  ] B=%-4ld N=%-4ld iters=%-7ld %8.4f ms/call  "
            "%8.3f us/event  %8.4f us/(event,var)\n",
            (long)B, (long)N_var, (long)n, ms_b, us_e, us_var);
    }

    // Combined.
    {
        auto [n, t] = run_loop(duration, [&]() {
            session.Run(Ort::RunOptions{nullptr}, binding);
            ssr::evaluate_joint_batched(
                spec, coefs_buf.data(),
                u_vars.data(), s_vars.data(),
                (int)B, (int)N_var, d_out.data());
        });
        double ms_b = 1000.0 * t / (double)n;
        double us_e = 1e6 * t / (double)n / (double)B;
        std::printf(
            "[trunk+poly ] B=%-4ld N=%-4ld iters=%-7ld %8.4f ms/call  "
            "%8.3f us/event\n",
            (long)B, (long)N_var, (long)n, ms_b, us_e);
    }

    return 0;
} catch (const std::exception& e) {
    std::fprintf(stderr, "exception: %s\n", e.what());
    return 1;
}
