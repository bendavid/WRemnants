// Single-core libtorch C++ benchmark for the shift+smear-reweight
// polyhead AOTI .pt2 produced by ``shift_smear_reweight_export.py``.
//
// Times the realistic deployment pattern: per event, *one* trunk
// forward at (y, c) -> joint_coefs, followed by N polynomial
// evaluations at (u, σ) values that share the same coefs. The trunk
// dominates per-event NN cost; the polynomial sum is cheap and
// amortizes across many shift / smear variations.
//
// Inputs:
//     <package.pt2> <indices.json>  [--batch B]  [--n-var N]
//                                   [--duration S]  [--cpu-pin C]
//
// Wrapper signature (matches ReweightPolyheadInference):
//     y_raw [B, n_features=3]  float32
//     c_raw [B, n_cond=5]      float32
// Output:
//     joint_coefs [B, n_basis]  float32
//
// Polynomial side:
//   per event b, evaluate d[b, k] for k = 0..N-1 at (u, σ) drawn
//   uniformly in [-1, 1] in standardized-target σ_y units. Reports
//   timings split into (trunk, poly) and combined.
//
// Build (sets RPATH so libtorch is picked up at runtime):
//
//   g++ -O2 -std=c++17 \
//       -I"$LIBTORCH/include" -I"$LIBTORCH/include/torch/csrc/api/include" \
//       -I scripts/corrections/muon_calibration \
//       -L"$LIBTORCH/lib" -Wl,-rpath,"$LIBTORCH/lib" \
//       bench_aoti_shift_smear_reweight.cpp \
//       -ltorch -ltorch_cpu -lc10 -lpthread \
//       -o bench_aoti_shift_smear_reweight

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
            "usage: %s <package.pt2> <indices.json> "
            "[--batch B] [--n-var N] [--duration S] [--cpu-pin C]\n",
            argv[0]);
        return 2;
    }
    const char* pt2_path  = argv[1];
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
    torch::set_num_threads(1);
    torch::set_num_interop_threads(1);

    // Load the polynomial sidecar (cheap, once).
    auto spec = ssr::load_basis_spec(json_path);

    // Load AOTI runner.
    torch::inductor::AOTIModelPackageLoader loader(pt2_path);
    std::fprintf(stderr,
        "loaded AOTI package: %s\n  B=%ld  N_var=%ld  n_basis=%d  "
        "n_features=%d  basis=%s\n",
        pt2_path, (long)B, (long)N_var, spec.n_basis,
        spec.n_features, spec.basis.c_str());

    const int64_t n_y = spec.n_features;
    const int64_t n_c = 5;  // ReweightPolyheadInference takes c[B, 5]
    auto opts = torch::TensorOptions()
                    .dtype(torch::kFloat32)
                    .device(torch::kCPU);
    auto y = torch::randn({B, n_y}, opts);
    auto c = torch::randn({B, n_c}, opts);
    std::vector<torch::Tensor> inputs{y, c};

    // Pre-allocate u, σ, and d buffers for the polynomial side (host-
    // owned, reused across iters; avoids heap allocation in the hot
    // loop). u and σ are filled once with random values in [-0.7, 0.7]
    // — the polynomial evaluation is data-independent in cost so the
    // exact distribution doesn't affect timings.
    std::mt19937 rng(0);
    std::uniform_real_distribution<float> uni(-0.7f, 0.7f);
    std::vector<float> u_vars  ((std::size_t)B * N_var * n_y);
    std::vector<float> s_vars  ((std::size_t)B * N_var * n_y);
    std::vector<float> d_out   ((std::size_t)B * N_var);
    for (auto& v : u_vars) v = uni(rng);
    for (auto& v : s_vars) v = uni(rng);

    // Warmup.
    for (int i = 0; i < 50; ++i) {
        auto out = loader.run(inputs);
        const auto& cf = out.at(0);
        ssr::evaluate_joint_batched(
            spec, cf.data_ptr<float>(),
            u_vars.data(), s_vars.data(),
            (int)B, (int)N_var, d_out.data());
    }

    // Three measurements: trunk-only, poly-only, combined.
    const int64_t min_iters = 5;
    auto run_loop = [&](double dur, auto&& body) {
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

    // Trunk only (one forward per call → coefs).
    {
        // Save coefs out of band so we can reuse for poly-only timing.
        torch::Tensor coefs_persist;
        auto [n, t] = run_loop(duration, [&]() {
            auto out = loader.run(inputs);
            coefs_persist = out.at(0);
        });
        double ms_b   = 1000.0 * t / (double)n;
        double us_e   = 1e6   * t / (double)n / (double)B;
        std::printf(
            "[trunk-only ] B=%-4ld iters=%-7ld %8.4f ms/call  "
            "%8.3f us/event  (%ld events)\n",
            (long)B, (long)n, ms_b, us_e, (long)B);
    }

    // Poly only (use a held coefs from one trunk call; pure C++ loop).
    {
        auto out0 = loader.run(inputs);
        auto coefs_t = out0.at(0).contiguous();
        const float* coefs_ptr = coefs_t.data_ptr<float>();
        auto [n, t] = run_loop(duration, [&]() {
            ssr::evaluate_joint_batched(
                spec, coefs_ptr,
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
            auto out = loader.run(inputs);
            const auto& cf = out.at(0);
            ssr::evaluate_joint_batched(
                spec, cf.data_ptr<float>(),
                u_vars.data(), s_vars.data(),
                (int)B, (int)N_var, d_out.data());
        });
        double ms_b   = 1000.0 * t / (double)n;
        double us_e   = 1e6   * t / (double)n / (double)B;
        std::printf(
            "[trunk+poly ] B=%-4ld N=%-4ld iters=%-7ld %8.4f ms/call  "
            "%8.3f us/event\n",
            (long)B, (long)N_var, (long)n, ms_b, us_e);
    }

    return 0;
} catch (const c10::Error& e) {
    std::fprintf(stderr, "torch error: %s\n", e.what());
    return 1;
} catch (const std::exception& e) {
    std::fprintf(stderr, "exception: %s\n", e.what());
    return 1;
}
