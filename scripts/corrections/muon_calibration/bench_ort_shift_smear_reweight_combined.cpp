// Single-core ONNX Runtime C++ benchmark for the *combined* shift+
// smear-reweight polyhead — the variant produced by
// ``shift_smear_reweight_export.py --combined-onnx-output ...`` that
// folds the polynomial evaluation into the model graph.
//
// Inference signature:
//     y_raw [B, n_features=3]    float32
//     c_raw [B, n_cond=5]        float32
//     u     [B, N_var, F]        float32
//     sigma [B, N_var, F]        float32
// Output:
//     d     [B, N_var]           float32
//
// The trunk runs once per event; the polynomial is evaluated at all
// N_var (u, σ) values in-graph.  No external evaluate_joint.hpp /
// JSON sidecar is needed at runtime.
//
// Build:
//   ORT_LIB=/opt/venv/lib/python3.13/site-packages/onnxruntime/capi
//   g++ -O2 -std=c++17 \
//       -I /usr/include/onnxruntime \
//       -L "$ORT_LIB" -Wl,-rpath,"$ORT_LIB" \
//       bench_ort_shift_smear_reweight_combined.cpp \
//       -lonnxruntime -lpthread \
//       -o bench_ort_ssr_combined

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
    if (argc < 2) {
        std::fprintf(stderr,
            "usage: %s <combined.onnx> "
            "[--batch B] [--n-var N] [--duration S] [--cpu-pin C] "
            "[--n-features F] [--n-cond NC]\n", argv[0]);
        return 2;
    }
    const char* onnx_path = argv[1];
    int64_t B     = 1;
    int64_t N_var = 32;
    double duration = 2.0;
    int cpu_pin = 0;
    int64_t n_y = 3;
    int64_t n_c = 5;
    for (int i = 2; i < argc; ++i) {
        if (std::strcmp(argv[i], "--batch") == 0 && i + 1 < argc) {
            B = std::atoll(argv[++i]);
        } else if (std::strcmp(argv[i], "--n-var") == 0 && i + 1 < argc) {
            N_var = std::atoll(argv[++i]);
        } else if (std::strcmp(argv[i], "--duration") == 0 && i + 1 < argc) {
            duration = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "--cpu-pin") == 0 && i + 1 < argc) {
            cpu_pin = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--n-features") == 0 && i + 1 < argc) {
            n_y = std::atoll(argv[++i]);
        } else if (std::strcmp(argv[i], "--n-cond") == 0 && i + 1 < argc) {
            n_c = std::atoll(argv[++i]);
        }
    }

    setenv("OMP_NUM_THREADS",      "1", 1);
    setenv("MKL_NUM_THREADS",      "1", 1);
    setenv("OPENBLAS_NUM_THREADS", "1", 1);
    pin_to_cpu(cpu_pin);

    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "bench_ssr_combined");
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
    if (n_inputs != 4 || n_outputs < 1) {
        std::fprintf(stderr,
            "expected 4 inputs / >=1 output; got %zu/%zu\n",
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

    std::fprintf(stderr,
        "model: %s\n"
        "inputs:  %s, %s, %s, %s\n"
        "B=%ld  N_var=%ld  n_features=%ld  n_cond=%ld\n",
        onnx_path,
        in_names[0].c_str(), in_names[1].c_str(),
        in_names[2].c_str(), in_names[3].c_str(),
        (long)B, (long)N_var, (long)n_y, (long)n_c);

    auto mem_info = Ort::MemoryInfo::CreateCpu(
        OrtArenaAllocator, OrtMemTypeDefault);

    std::mt19937 rng(0);
    std::normal_distribution<float> gauss(0.0f, 1.0f);
    std::uniform_real_distribution<float> uni(-0.7f, 0.7f);

    std::vector<float> y_buf((std::size_t)B * n_y);
    std::vector<float> c_buf((std::size_t)B * n_c);
    std::vector<float> u_buf((std::size_t)B * N_var * n_y);
    std::vector<float> s_buf((std::size_t)B * N_var * n_y);
    std::vector<float> d_buf((std::size_t)B * N_var);
    for (auto& v : y_buf) v = gauss(rng);
    for (auto& v : c_buf) v = gauss(rng);
    for (auto& v : u_buf) v = uni(rng);
    for (auto& v : s_buf) v = uni(rng);

    std::array<int64_t, 2> y_shape{B, n_y};
    std::array<int64_t, 2> c_shape{B, n_c};
    std::array<int64_t, 3> u_shape{B, N_var, n_y};
    std::array<int64_t, 3> s_shape{B, N_var, n_y};
    std::array<int64_t, 2> d_shape{B, N_var};

    Ort::Value y_t = Ort::Value::CreateTensor<float>(
        mem_info, y_buf.data(), y_buf.size(),
        y_shape.data(), y_shape.size());
    Ort::Value c_t = Ort::Value::CreateTensor<float>(
        mem_info, c_buf.data(), c_buf.size(),
        c_shape.data(), c_shape.size());
    Ort::Value u_t = Ort::Value::CreateTensor<float>(
        mem_info, u_buf.data(), u_buf.size(),
        u_shape.data(), u_shape.size());
    Ort::Value s_t = Ort::Value::CreateTensor<float>(
        mem_info, s_buf.data(), s_buf.size(),
        s_shape.data(), s_shape.size());
    Ort::Value d_t = Ort::Value::CreateTensor<float>(
        mem_info, d_buf.data(), d_buf.size(),
        d_shape.data(), d_shape.size());

    Ort::IoBinding binding(session);
    binding.BindInput (in_names[0].c_str(), y_t);
    binding.BindInput (in_names[1].c_str(), c_t);
    binding.BindInput (in_names[2].c_str(), u_t);
    binding.BindInput (in_names[3].c_str(), s_t);
    binding.BindOutput(out_names[0].c_str(), d_t);

    // Warmup.
    for (int i = 0; i < 50; ++i) {
        session.Run(Ort::RunOptions{nullptr}, binding);
    }

    const int64_t min_iters = 5;
    int64_t n_iters = 0;
    auto t0 = std::chrono::steady_clock::now();
    auto deadline = t0 + std::chrono::duration_cast<
        std::chrono::steady_clock::duration>(
        std::chrono::duration<double>(duration));
    while (true) {
        session.Run(Ort::RunOptions{nullptr}, binding);
        ++n_iters;
        if (n_iters >= min_iters &&
            std::chrono::steady_clock::now() >= deadline) break;
    }
    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    double ms_b = 1000.0 * elapsed / (double)n_iters;
    double us_e = 1e6 * elapsed / (double)n_iters / (double)B;
    double us_v = us_e / (double)N_var;
    std::printf(
        "[combined ] B=%-4ld N=%-4ld iters=%-7ld %8.4f ms/call  "
        "%8.3f us/event  %8.4f us/(event,var)\n",
        (long)B, (long)N_var, (long)n_iters, ms_b, us_e, us_v);
    return 0;
} catch (const std::exception& e) {
    std::fprintf(stderr, "exception: %s\n", e.what());
    return 1;
}
