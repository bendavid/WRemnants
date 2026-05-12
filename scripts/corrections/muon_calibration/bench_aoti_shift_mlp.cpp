// Single-core libtorch C++ benchmark for the shift-reweight-MLP AOTI .pt2
// produced by ``train_shift_reweight_mlp_export.py``.
//
// Each .pt2 was compiled for a specific static batch size; this binary
// takes the path + batch size and times forward calls. Pin the process
// to a single CPU core, disable threading, run a warmup, then measure
// for ``--duration`` seconds.
//
// Inputs:    <package.pt2> <batch> [--duration S] [--cpu-pin C]
//
// Wrapper signature (matches ShiftReweightInference):
//     y_raw [B, 3]  float32     (raw target, e.g. r_kappa,dlambda,dphi)
//     c_raw [B, 5]  float32     (raw conditioning)
//     u_raw [B, 3]  float32     (y-space shift in physical units)
// Output:
//     log_r [B]     float32     (= log p(y - u, c) / p(y, c))
//
// Example build (adjust paths to your libtorch install):
//     g++ -O2 -std=c++17 \
//         -I"$LIBTORCH/include" -I"$LIBTORCH/include/torch/csrc/api/include" \
//         -L"$LIBTORCH/lib" \
//         -Wl,-rpath,"$LIBTORCH/lib" \
//         bench_aoti_shift_mlp.cpp \
//         -ltorch -ltorch_cpu -lc10 -lpthread \
//         -o bench_aoti_shift_mlp

#define _GNU_SOURCE
#include <sched.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <torch/torch.h>
#include <torch/csrc/inductor/aoti_package/model_package_loader.h>

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
            "usage: %s <package.pt2> <batch> [--duration S] [--cpu-pin C]\n",
            argv[0]);
        return 2;
    }
    const char* pt2_path = argv[1];
    int64_t B = std::atoll(argv[2]);
    double duration = 2.0;
    int cpu_pin = 0;
    for (int i = 3; i < argc; ++i) {
        if (std::strcmp(argv[i], "--duration") == 0 && i + 1 < argc) {
            duration = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "--cpu-pin") == 0 && i + 1 < argc) {
            cpu_pin = std::atoi(argv[++i]);
        }
    }

    setenv("OMP_NUM_THREADS",     "1", 1);
    setenv("MKL_NUM_THREADS",     "1", 1);
    setenv("OPENBLAS_NUM_THREADS", "1", 1);
    pin_to_cpu(cpu_pin);
    torch::set_num_threads(1);
    torch::set_num_interop_threads(1);

    torch::inductor::AOTIModelPackageLoader loader(pt2_path);
    std::fprintf(stderr, "loaded AOTI package: %s  (B=%ld)\n",
                 pt2_path, (long)B);

    const int64_t n_y = 3, n_c = 5, n_u = 3;
    auto opts = torch::TensorOptions()
                    .dtype(torch::kFloat32)
                    .device(torch::kCPU);
    auto y = torch::randn({B, n_y}, opts);
    auto c = torch::randn({B, n_c}, opts);
    auto u = torch::randn({B, n_u}, opts) * 0.1f;
    std::vector<torch::Tensor> inputs{y, c, u};

    // Warmup.
    for (int i = 0; i < 50; ++i) {
        auto out = loader.run(inputs);
        (void)out;
    }

    const int64_t min_iters = 5;
    int64_t n_iters = 0;
    auto t0 = std::chrono::steady_clock::now();
    auto deadline = t0 + std::chrono::duration_cast<
        std::chrono::steady_clock::duration>(
        std::chrono::duration<double>(duration));
    while (true) {
        auto out = loader.run(inputs);
        (void)out;
        ++n_iters;
        if (n_iters >= min_iters &&
            std::chrono::steady_clock::now() >= deadline) break;
    }
    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    double ms_b = 1000.0 * elapsed / (double)n_iters;
    double us_e = 1e6 * elapsed / (double)n_iters / (double)B;
    double ev_s = (double)n_iters * (double)B / elapsed;
    std::printf("%8s %10s %10s %10s %14s\n",
                "batch", "iters", "ms/batch", "us/event", "events/s");
    std::printf("%8ld %10ld %10.4f %10.3f %14.1f\n",
                (long)B, (long)n_iters, ms_b, us_e, ev_s);
    return 0;
} catch (const c10::Error& e) {
    std::fprintf(stderr, "torch error: %s\n", e.what());
    return 1;
} catch (const std::exception& e) {
    std::fprintf(stderr, "exception: %s\n", e.what());
    return 1;
}
