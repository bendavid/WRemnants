// Single-core libtorch C++ benchmark for the *combined* shift+smear-
// reweight polyhead AOTI .pt2 produced by
// ``shift_smear_reweight_export.py --combined-output ...``.
//
// Inference signature (matches CombinedInference):
//     y_raw [B, n_features=3]     float32
//     c_raw [B, n_cond=5]         float32
//     u     [B, N_var, F]         float32
//     sigma [B, N_var, F]         float32
// Output:
//     d     [B, N_var]            float32
//
// The trunk runs once per event; the polynomial sum is evaluated at
// all N_var (u, σ) values in-graph. No external evaluate_joint.hpp
// needed at runtime.
//
// Build:
//   LIBTORCH=/opt/venv/lib/python3.13/site-packages/torch
//   g++ -O2 -std=c++17 \
//       -I"$LIBTORCH/include" -I"$LIBTORCH/include/torch/csrc/api/include" \
//       -L"$LIBTORCH/lib" -Wl,-rpath,"$LIBTORCH/lib" \
//       bench_aoti_shift_smear_reweight_combined.cpp \
//       -ltorch -ltorch_cpu -lc10 -lpthread \
//       -o bench_aoti_ssr_combined

#define _GNU_SOURCE
#include <sched.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
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
    if (argc < 2) {
        std::fprintf(stderr,
            "usage: %s <combined.pt2> "
            "[--batch B] [--n-var N] [--duration S] [--cpu-pin C] "
            "[--n-features F] [--n-cond NC]\n", argv[0]);
        return 2;
    }
    const char* pt2_path = argv[1];
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
    torch::set_num_threads(1);
    torch::set_num_interop_threads(1);

    torch::inductor::AOTIModelPackageLoader loader(pt2_path);
    std::fprintf(stderr,
        "loaded combined AOTI: %s\n  B=%ld  N_var=%ld  "
        "n_features=%ld  n_cond=%ld\n",
        pt2_path, (long)B, (long)N_var, (long)n_y, (long)n_c);

    auto opts = torch::TensorOptions()
                    .dtype(torch::kFloat32)
                    .device(torch::kCPU);
    auto y = torch::randn({B, n_y}, opts);
    auto c = torch::randn({B, n_c}, opts);
    auto u = torch::empty({B, N_var, n_y}, opts).uniform_(-0.7f, 0.7f);
    auto s = torch::empty({B, N_var, n_y}, opts).uniform_(-0.7f, 0.7f);
    std::vector<torch::Tensor> inputs{y, c, u, s};

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
    double us_v = us_e / (double)N_var;
    std::printf(
        "[combined ] B=%-4ld N=%-4ld iters=%-7ld %8.4f ms/call  "
        "%8.3f us/event  %8.4f us/(event,var)\n",
        (long)B, (long)N_var, (long)n_iters, ms_b, us_e, us_v);
    return 0;
} catch (const c10::Error& e) {
    std::fprintf(stderr, "torch error: %s\n", e.what());
    return 1;
} catch (const std::exception& e) {
    std::fprintf(stderr, "exception: %s\n", e.what());
    return 1;
}
