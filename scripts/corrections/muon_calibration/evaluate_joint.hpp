// Header-only polynomial evaluator for the shift+smear-reweight
// polyhead's joint basis. Pairs with the JSON sidecar emitted by
// shift_smear_reweight_export.py:
//
//     {
//       "n_features":  3,
//       "n_basis":     122,
//       "alpha_degs":  [[k_u_0, k_u_1, k_u_2], ...],   // [n_basis, F]
//       "beta_degs":   [[k_s_0, k_s_1, k_s_2], ...],
//       "basis":       "monomial" | "chebyshev",
//       "basis_scale_u":     <float>,
//       "basis_scale_sigma": <float>,
//       "max_deg_u":   <int>,
//       "max_deg_sigma": <int>,
//       ...
//     }
//
// The C++ benches only need n_features, n_basis, alpha_degs,
// beta_degs, basis, basis_scale_u, basis_scale_sigma, max_deg_u,
// max_deg_sigma — the parser is therefore intentionally minimal and
// hand-rolled (no third-party JSON dep on the bench harness side).
//
// Evaluation:
//
//     d[b, k] = Σ_i coefs[b, i] · phi_i(u[b, k, :], σ[b, k, :])
//
// where phi_i = Π_j basis(u_j; α_{i,j}) · basis(σ_j; β_{i,j}).
// "monomial" uses x^n; "chebyshev" uses tensor-product T_n(x / scale)
// with the structural zero-anchor (subtract Π T_n(0)) so r=1 at
// (u, σ)=(0, 0).
//
// All call signatures take strided float32 arrays in row-major
// layout (C contiguous):
//
//     coefs    [B, n_basis]
//     u        [B, N_var, n_features]
//     sigma    [B, N_var, n_features]
//     d_out    [B, N_var]
//
// No dynamic allocation in the hot path beyond a small per-event
// scratch buffer for the per-axis power tables.

#pragma once

#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ssr {

struct BasisSpec {
    int n_features = 0;
    int n_basis    = 0;
    int max_deg_u  = 0;
    int max_deg_sigma = 0;
    std::string basis = "monomial";
    float scale_u     = 1.0f;
    float scale_sigma = 1.0f;
    // Row-major [n_basis, n_features].
    std::vector<int> alpha_degs;
    std::vector<int> beta_degs;
    // Per-side Chebyshev origin constants for the doubly-centered
    // factorization phi[k] = (Π T_α(u_norm) − cheb_u_const[k])
    //                     · (Π T_β(σ_norm) − cheb_v_const[k]).
    // cheb_u_const[k] = Π T_α(0)  if α nonempty else 0;
    // cheb_v_const[k] = Π T_β(0)  if β nonempty else 0.
    // Per-side anchoring guarantees pure-u / pure-σ / cross terms
    // vanish at u=0 / σ=0 / either, matching the structural priors
    // and the per-mode trainer losses. For monomial these are
    // unused.
    std::vector<float> cheb_u_const;
    std::vector<float> cheb_v_const;
};

// ---------------------------------------------------------------------------
// Minimal JSON sidecar parser.
//
// Schema is fixed and hand-controlled: top-level object with the
// keys listed above. We only support the shapes we emit.
// ---------------------------------------------------------------------------
namespace detail {

inline void _skip_ws(const std::string& s, std::size_t& p) {
    while (p < s.size() &&
           (s[p] == ' ' || s[p] == '\t' || s[p] == '\n' || s[p] == '\r')) {
        ++p;
    }
}

inline std::string _parse_string(const std::string& s, std::size_t& p) {
    _skip_ws(s, p);
    if (p >= s.size() || s[p] != '"') {
        throw std::runtime_error("JSON: expected '\"' at offset " +
                                 std::to_string(p));
    }
    ++p;
    std::string out;
    while (p < s.size() && s[p] != '"') {
        if (s[p] == '\\' && p + 1 < s.size()) {
            out.push_back(s[p + 1]);
            p += 2;
        } else {
            out.push_back(s[p]);
            ++p;
        }
    }
    if (p >= s.size()) throw std::runtime_error("JSON: unterminated string");
    ++p;  // closing quote
    return out;
}

inline double _parse_number(const std::string& s, std::size_t& p) {
    _skip_ws(s, p);
    std::size_t start = p;
    if (p < s.size() && (s[p] == '-' || s[p] == '+')) ++p;
    while (p < s.size() &&
           ((s[p] >= '0' && s[p] <= '9') || s[p] == '.' ||
            s[p] == 'e' || s[p] == 'E' || s[p] == '+' || s[p] == '-')) {
        ++p;
    }
    return std::stod(s.substr(start, p - start));
}

inline void _expect(const std::string& s, std::size_t& p, char c) {
    _skip_ws(s, p);
    if (p >= s.size() || s[p] != c) {
        throw std::runtime_error(std::string("JSON: expected '") + c +
                                 "' at offset " + std::to_string(p));
    }
    ++p;
}

inline std::vector<int> _parse_int_array(const std::string& s,
                                         std::size_t& p) {
    std::vector<int> out;
    _expect(s, p, '[');
    _skip_ws(s, p);
    if (p < s.size() && s[p] == ']') {
        ++p;
        return out;
    }
    while (true) {
        out.push_back(static_cast<int>(_parse_number(s, p)));
        _skip_ws(s, p);
        if (p < s.size() && s[p] == ',') { ++p; continue; }
        break;
    }
    _expect(s, p, ']');
    return out;
}

inline std::vector<std::vector<int>>
_parse_int_array_2d(const std::string& s, std::size_t& p) {
    std::vector<std::vector<int>> out;
    _expect(s, p, '[');
    _skip_ws(s, p);
    if (p < s.size() && s[p] == ']') {
        ++p;
        return out;
    }
    while (true) {
        out.push_back(_parse_int_array(s, p));
        _skip_ws(s, p);
        if (p < s.size() && s[p] == ',') { ++p; continue; }
        break;
    }
    _expect(s, p, ']');
    return out;
}

// Skip past a JSON value that we don't care about (so the parser is
// resilient to extra keys like target_mean / cond_std that the bench
// doesn't need).
inline void _skip_value(const std::string& s, std::size_t& p) {
    _skip_ws(s, p);
    if (p >= s.size()) throw std::runtime_error("JSON: unexpected EOF");
    char c = s[p];
    if (c == '"') { (void)_parse_string(s, p); return; }
    if (c == '[') {
        ++p;
        int depth = 1;
        bool in_str = false;
        while (p < s.size() && depth > 0) {
            if (in_str) {
                if (s[p] == '\\' && p + 1 < s.size()) { p += 2; continue; }
                if (s[p] == '"') in_str = false;
            } else {
                if (s[p] == '"') in_str = true;
                else if (s[p] == '[') ++depth;
                else if (s[p] == ']') --depth;
            }
            ++p;
        }
        return;
    }
    if (c == '{') {
        ++p;
        int depth = 1;
        bool in_str = false;
        while (p < s.size() && depth > 0) {
            if (in_str) {
                if (s[p] == '\\' && p + 1 < s.size()) { p += 2; continue; }
                if (s[p] == '"') in_str = false;
            } else {
                if (s[p] == '"') in_str = true;
                else if (s[p] == '{') ++depth;
                else if (s[p] == '}') --depth;
            }
            ++p;
        }
        return;
    }
    // number / bool / null
    while (p < s.size() && s[p] != ',' && s[p] != '}' && s[p] != ']') {
        ++p;
    }
}

}  // namespace detail

inline BasisSpec load_basis_spec(const std::string& json_path) {
    std::ifstream f(json_path);
    if (!f) {
        throw std::runtime_error("cannot open " + json_path);
    }
    std::stringstream ss;
    ss << f.rdbuf();
    std::string s = ss.str();
    std::size_t p = 0;

    BasisSpec spec;
    detail::_expect(s, p, '{');
    detail::_skip_ws(s, p);
    std::vector<std::vector<int>> aa, bb;
    while (true) {
        detail::_skip_ws(s, p);
        if (p < s.size() && s[p] == '}') { ++p; break; }
        std::string key = detail::_parse_string(s, p);
        detail::_expect(s, p, ':');
        if (key == "n_features") {
            spec.n_features = static_cast<int>(detail::_parse_number(s, p));
        } else if (key == "n_basis") {
            spec.n_basis = static_cast<int>(detail::_parse_number(s, p));
        } else if (key == "max_deg_u") {
            spec.max_deg_u = static_cast<int>(detail::_parse_number(s, p));
        } else if (key == "max_deg_sigma") {
            spec.max_deg_sigma =
                static_cast<int>(detail::_parse_number(s, p));
        } else if (key == "basis") {
            spec.basis = detail::_parse_string(s, p);
        } else if (key == "basis_scale_u") {
            spec.scale_u = static_cast<float>(detail::_parse_number(s, p));
        } else if (key == "basis_scale_sigma") {
            spec.scale_sigma =
                static_cast<float>(detail::_parse_number(s, p));
        } else if (key == "alpha_degs") {
            aa = detail::_parse_int_array_2d(s, p);
        } else if (key == "beta_degs") {
            bb = detail::_parse_int_array_2d(s, p);
        } else {
            detail::_skip_value(s, p);
        }
        detail::_skip_ws(s, p);
        if (p < s.size() && s[p] == ',') { ++p; continue; }
    }
    if ((int)aa.size() != spec.n_basis ||
        (int)bb.size() != spec.n_basis) {
        throw std::runtime_error(
            "JSON: alpha/beta_degs row count != n_basis");
    }
    spec.alpha_degs.resize((std::size_t)spec.n_basis * spec.n_features);
    spec.beta_degs.resize((std::size_t)spec.n_basis * spec.n_features);
    for (int i = 0; i < spec.n_basis; ++i) {
        if ((int)aa[i].size() != spec.n_features ||
            (int)bb[i].size() != spec.n_features) {
            throw std::runtime_error(
                "JSON: alpha/beta_degs row width != n_features");
        }
        for (int j = 0; j < spec.n_features; ++j) {
            spec.alpha_degs[i * spec.n_features + j] = aa[i][j];
            spec.beta_degs [i * spec.n_features + j] = bb[i][j];
        }
    }
    spec.cheb_u_const.assign((std::size_t)spec.n_basis, 0.0f);
    spec.cheb_v_const.assign((std::size_t)spec.n_basis, 0.0f);
    if (spec.basis == "chebyshev") {
        // T_n(0) = 0 if n is odd, (-1)^(n/2) if n is even, 1 if n == 0.
        auto T_at_zero = [](int n) -> float {
            if (n == 0) return 1.0f;
            if (n % 2 == 1) return 0.0f;
            return (((n / 2) % 2) == 0) ? 1.0f : -1.0f;
        };
        for (int i = 0; i < spec.n_basis; ++i) {
            bool a_nonempty = false, b_nonempty = false;
            for (int j = 0; j < spec.n_features; ++j) {
                if (spec.alpha_degs[i * spec.n_features + j] > 0)
                    a_nonempty = true;
                if (spec.beta_degs [i * spec.n_features + j] > 0)
                    b_nonempty = true;
            }
            if (a_nonempty) {
                float cu = 1.0f;
                for (int j = 0; j < spec.n_features; ++j) {
                    int au = spec.alpha_degs[i * spec.n_features + j];
                    if (au > 0) cu *= T_at_zero(au);
                    if (cu == 0.0f) break;
                }
                spec.cheb_u_const[i] = cu;
            }
            if (b_nonempty) {
                float cv = 1.0f;
                for (int j = 0; j < spec.n_features; ++j) {
                    int bu = spec.beta_degs [i * spec.n_features + j];
                    if (bu > 0) cv *= T_at_zero(bu);
                    if (cv == 0.0f) break;
                }
                spec.cheb_v_const[i] = cv;
            }
        }
    }
    return spec;
}

// ---------------------------------------------------------------------------
// Polynomial evaluation.
//
// Per-event power tables:
//   pu[F, max_deg_u + 1]      : pu[j, k] = u_j^k or T_k(u_j / scale_u)
//   ps[F, max_deg_sigma + 1]  : ps[j, k] = σ_j^k or T_k(σ_j / scale_sigma)
//
// Then phi_i = Π_j pu[j, α_{i,j}] · ps[j, β_{i,j}], possibly minus
// the Chebyshev origin constant.
// ---------------------------------------------------------------------------

// One-event evaluation. ``coefs`` shape [n_basis], ``u``/``sigma``
// shape [N_var, F]. ``d_out`` filled with N_var values.
inline void evaluate_joint_event(
    const BasisSpec& spec,
    const float* __restrict__ coefs,    // [n_basis]
    const float* __restrict__ u_vars,   // [N_var, F]
    const float* __restrict__ sigma_vars, // [N_var, F]
    int N_var,
    float* __restrict__ d_out,          // [N_var]
    // Reusable scratch — pass nullptr to allocate inline.
    float* __restrict__ pow_u_buf = nullptr,
    float* __restrict__ pow_s_buf = nullptr
) {
    const int F  = spec.n_features;
    const int Nb = spec.n_basis;
    const int Du = spec.max_deg_u;
    const int Ds = spec.max_deg_sigma;

    std::vector<float> pow_u_owned, pow_s_owned;
    if (!pow_u_buf) {
        pow_u_owned.assign((std::size_t)F * (Du + 1), 0.0f);
        pow_u_buf = pow_u_owned.data();
    }
    if (!pow_s_buf) {
        pow_s_owned.assign((std::size_t)F * (Ds + 1), 0.0f);
        pow_s_buf = pow_s_owned.data();
    }

    const bool is_cheb = (spec.basis == "chebyshev");
    const float inv_scale_u = is_cheb ? (1.0f / spec.scale_u) : 1.0f;
    const float inv_scale_s = is_cheb ? (1.0f / spec.scale_sigma) : 1.0f;

    for (int k = 0; k < N_var; ++k) {
        const float* u_k = u_vars     + (std::size_t)k * F;
        const float* s_k = sigma_vars + (std::size_t)k * F;

        // Build per-axis power / Chebyshev tables.
        for (int j = 0; j < F; ++j) {
            float* pu = pow_u_buf + (std::size_t)j * (Du + 1);
            float* ps = pow_s_buf + (std::size_t)j * (Ds + 1);
            const float xu = u_k[j] * inv_scale_u;
            const float xs = s_k[j] * inv_scale_s;
            pu[0] = 1.0f;
            ps[0] = 1.0f;
            if (Du >= 1) pu[1] = is_cheb ? xu : xu;
            if (Ds >= 1) ps[1] = is_cheb ? xs : xs;
            if (is_cheb) {
                for (int n = 2; n <= Du; ++n) {
                    pu[n] = 2.0f * xu * pu[n - 1] - pu[n - 2];
                }
                for (int n = 2; n <= Ds; ++n) {
                    ps[n] = 2.0f * xs * ps[n - 1] - ps[n - 2];
                }
            } else {
                for (int n = 2; n <= Du; ++n) {
                    pu[n] = pu[n - 1] * xu;
                }
                for (int n = 2; n <= Ds; ++n) {
                    ps[n] = ps[n - 1] * xs;
                }
            }
        }

        // Accumulate d = Σ_i coefs[i] · phi_i.
        // For chebyshev we use the doubly-centered factorization
        //   phi = (U_full - cheb_u_const[i]) * (V_full - cheb_v_const[i])
        // where U_full = Π pu[j, α_j], V_full = Π ps[j, β_j], and the
        // per-side constants are zero on whichever side is empty so
        // pure-u / pure-σ multi-indices keep their singly-centered
        // forms (V-0)=1 or (U-0)=1.
        float d = 0.0f;
        for (int i = 0; i < Nb; ++i) {
            const int* a = spec.alpha_degs.data() + (std::size_t)i * F;
            const int* b = spec.beta_degs .data() + (std::size_t)i * F;
            float u_factor = 1.0f;
            float v_factor = 1.0f;
            for (int j = 0; j < F; ++j) {
                if (a[j]) u_factor *= pow_u_buf[(std::size_t)j * (Du + 1) + a[j]];
                if (b[j]) v_factor *= pow_s_buf[(std::size_t)j * (Ds + 1) + b[j]];
            }
            float phi;
            if (is_cheb) {
                phi = (u_factor - spec.cheb_u_const[i])
                    * (v_factor - spec.cheb_v_const[i]);
            } else {
                phi = u_factor * v_factor;
            }
            d += coefs[i] * phi;
        }
        d_out[k] = d;
    }
}

// Batched evaluation. Inputs row-major:
//   coefs   [B, n_basis]
//   u       [B, N_var, F]
//   sigma   [B, N_var, F]
//   d_out   [B, N_var]
inline void evaluate_joint_batched(
    const BasisSpec& spec,
    const float* coefs,
    const float* u_vars,
    const float* sigma_vars,
    int B, int N_var,
    float* d_out
) {
    const int F  = spec.n_features;
    const int Nb = spec.n_basis;
    const int Du = spec.max_deg_u;
    const int Ds = spec.max_deg_sigma;

    // Per-call scratch reused across events.
    std::vector<float> pow_u((std::size_t)F * (Du + 1));
    std::vector<float> pow_s((std::size_t)F * (Ds + 1));

    for (int b = 0; b < B; ++b) {
        evaluate_joint_event(
            spec,
            coefs    + (std::size_t)b * Nb,
            u_vars   + (std::size_t)b * N_var * F,
            sigma_vars + (std::size_t)b * N_var * F,
            N_var,
            d_out    + (std::size_t)b * N_var,
            pow_u.data(), pow_s.data()
        );
    }
}

}  // namespace ssr
