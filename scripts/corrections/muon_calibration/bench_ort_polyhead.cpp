// Single-core ONNX Runtime C++ benchmark for the flow polyhead ONNX
// produced by ``flow_polyhead_export_onnx.py``.
//
// Builds a session, preallocates input/output tensors, binds them via
// IoBinding, and runs a warm loop for a fixed wall-clock duration at
// each batch size. Pins to a single CPU core; disables intra/inter-op
// threading to make the numbers comparable to the single-event narf
// runtime.
//
// Inputs:    <model.onnx> [--duration S] [--cpu-pin C] [--batches B1 B2 ...]
//
// Wrapper signature (matches FlowPolyheadInference):
//     y_raw [B, 3]  float32
//     c_raw [B, 5]  float32
// Output:
//     joint_coefs [B, n_basis] float32
//
// Example build (adjust ORT path):
//     g++ -O2 -std=c++17 \
//         -I"$ORT/include" -L"$ORT/lib" -Wl,-rpath,"$ORT/lib" \
//         bench_ort_polyhead.cpp \
//         -lonnxruntime -lpthread \
//         -o bench_ort_polyhead

#define _GNU_SOURCE
#include <sched.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
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
            "usage: %s <model.onnx> [--duration S] [--cpu-pin C] "
            "[--batches B1 B2 ...]\n", argv[0]);
        return 2;
    }
    const char* onnx_path = argv[1];
    double duration = 2.0;
    int cpu_pin = 0;
    std::vector<int64_t> batch_sizes = {
        1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
        1024, 2048, 4096, 8192, 16384,
    };
    bool batches_overridden = false;
    for (int i = 2; i < argc; ++i) {
        if (std::strcmp(argv[i], "--duration") == 0 && i + 1 < argc) {
            duration = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "--cpu-pin") == 0 && i + 1 < argc) {
            cpu_pin = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--batches") == 0) {
            if (!batches_overridden) batch_sizes.clear();
            batches_overridden = true;
            while (i + 1 < argc && argv[i + 1][0] != '-') {
                batch_sizes.push_back(std::atoll(argv[++i]));
            }
        }
    }

    setenv("OMP_NUM_THREADS",     "1", 1);
    setenv("MKL_NUM_THREADS",     "1", 1);
    setenv("OPENBLAS_NUM_THREADS", "1", 1);
    pin_to_cpu(cpu_pin);

    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "bench_polyhead");
    Ort::SessionOptions opts;
    opts.SetIntraOpNumThreads(1);
    opts.SetInterOpNumThreads(1);
    opts.SetExecutionMode(ORT_SEQUENTIAL);
    opts.SetGraphOptimizationLevel(
        GraphOptimizationLevel::ORT_ENABLE_ALL);
    opts.EnableCpuMemArena();

    Ort::Session session(env, onnx_path, opts);
    Ort::AllocatorWithDefaultOptions alloc;

    const size_t n_inputs  = session.GetInputCount();
    const size_t n_outputs = session.GetOutputCount();
    if (n_inputs != 2 || n_outputs < 1) {
        std::fprintf(stderr,
            "expected 2 inputs and >=1 output; got %zu / %zu\n",
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

    auto infer_dim1 = [&](size_t idx, int64_t fallback) -> int64_t {
        Ort::TypeInfo ti = session.GetInputTypeInfo(idx);
        auto info = ti.GetTensorTypeAndShapeInfo();
        const size_t rank = info.GetDimensionsCount();
        if (rank < 2) return fallback;
        std::vector<int64_t> shape(rank);
        info.GetDimensions(shape.data(), rank);
        return shape[1] > 0 ? shape[1] : fallback;
    };
    const int64_t n_y = infer_dim1(0, 3);
    const int64_t n_c = infer_dim1(1, 5);

    std::fprintf(stderr,
        "model: %s\ninputs: %s[N,%ld], %s[N,%ld]  outputs: %zu\n",
        onnx_path,
        in_names[0].c_str(), (long)n_y,
        in_names[1].c_str(), (long)n_c,
        n_outputs);

    auto mem_info = Ort::MemoryInfo::CreateCpu(
        OrtArenaAllocator, OrtMemTypeDefault);

    std::mt19937 rng(0);
    std::normal_distribution<float> gauss(0.0f, 1.0f);

    std::printf("%8s %10s %10s %10s %14s\n",
                "batch", "iters", "ms/batch", "us/event", "events/s");

    for (int64_t B : batch_sizes) {
        std::vector<float> y_buf(B * n_y);
        std::vector<float> c_buf(B * n_c);
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

        Ort::IoBinding binding(session);
        binding.BindInput(in_names[0].c_str(), y_tensor);
        binding.BindInput(in_names[1].c_str(), c_tensor);
        for (auto& oname : out_names) {
            binding.BindOutput(oname.c_str(), mem_info);
        }

        for (int i = 0; i < 10; ++i) {
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
        double ev_s = (double)n_iters * (double)B / elapsed;

        std::printf("%8ld %10ld %10.4f %10.3f %14.1f\n",
                    (long)B, (long)n_iters, ms_b, us_e, ev_s);
        std::fflush(stdout);
    }

    return 0;
} catch (const std::exception& e) {
    std::fprintf(stderr, "exception: %s\n", e.what());
    return 1;
} catch (...) {
    std::fprintf(stderr, "unknown exception\n");
    return 1;
}
