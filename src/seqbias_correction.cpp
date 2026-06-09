// Beta-Uniform mixture model for sequence bias correction
//
// Two exported functions (accessible as nomeR:::.<name>):
//
//   .fit_one_context_cpp(pos_data, neg_data, ...)
//     Joint constrained EM for one sequence context.
//     Returns alpha/beta/eps for pos and neg components plus convergence flags.
//
//   .correct_posterior_cpp(raw_p, param_idx, alpha_pos, beta_pos, eps_pos,
//                          alpha_neg, beta_neg, eps_neg, pi_pos)
//     Computes Bayesian posterior probabilities and applies pool-adjacent-
//     violators (PAV) isotonic regression per context. Single C++ pass.
//     Returns corrected probabilities (same length as raw_p).

#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
using namespace Rcpp;

// ---------------------------------------------------------------
// Scalar wrappers around R C-API special functions
// ---------------------------------------------------------------
static inline double digamma_c (double x) { return Rf_digamma(x);  }
static inline double trigamma_c(double x) { return Rf_trigamma(x); }
static inline double dbeta_c   (double x, double a, double b) {
    return R::dbeta(x, a, b, 0);
}

// ---------------------------------------------------------------
// Weighted sufficient statistics for Beta MLE
//   t1 = E_w[log x],   t2 = E_w[log(1-x)]
// Falls back to unweighted when sum(w) is negligible.
// ---------------------------------------------------------------
static void suff_stats(const NumericVector& lx,
                       const NumericVector& l1mx,
                       const NumericVector& w,
                       double& t1, double& t2) {
    int n = lx.size();
    double sw = 0, s1 = 0, s2 = 0;
    for (int i = 0; i < n; i++) { sw += w[i]; s1 += w[i]*lx[i]; s2 += w[i]*l1mx[i]; }
    if (sw < 1e-300) {
        s1 = s2 = 0;
        for (int i = 0; i < n; i++) { s1 += lx[i]; s2 += l1mx[i]; }
        t1 = s1/n; t2 = s2/n;
    } else { t1 = s1/sw; t2 = s2/sw; }
}

// ---------------------------------------------------------------
// Weighted Beta MLE via closed-form 2x2 Newton-Raphson
// Not exported to R; called from fit_beta_unif_em.
// ---------------------------------------------------------------
static NumericVector mle_beta_nr(const NumericVector& lx,
                                  const NumericVector& l1mx,
                                  const NumericVector& w,
                                  double a, double b,
                                  int max_iter = 30, double tol = 1e-9) {
    double t1, t2;
    suff_stats(lx, l1mx, w, t1, t2);
    for (int i = 0; i < max_iter; i++) {
        double ab  = a + b;
        double pab = digamma_c(ab);
        double g1  = pab - digamma_c(a) + t1;
        double g2  = pab - digamma_c(b) + t2;
        if (!std::isfinite(g1) || !std::isfinite(g2)) break;
        if (std::max(std::abs(g1), std::abs(g2)) < tol) break;
        double ta = trigamma_c(a), tb = trigamma_c(b), tab = trigamma_c(ab);
        double h11 = tab-ta, h12 = tab, h22 = tab-tb;
        double det = h11*h22 - h12*h12;
        if (!std::isfinite(det) || std::abs(det) < 1e-300) break;
        double an = a - (h22*g1 - h12*g2)/det;
        double bn = b - (-h12*g1 + h11*g2)/det;
        if (!std::isfinite(an) || !std::isfinite(bn) || an <= 0 || bn <= 0) break;
        if (std::max(std::abs(an-a), std::abs(bn-b)) < tol) { a=an; b=bn; break; }
        a = an; b = bn;
    }
    return NumericVector::create(Named("alpha")=a, Named("beta")=b);
}

// ---------------------------------------------------------------
// EM for a single Beta-Uniform mixture
//   f(x|a,b,e) = (1-e)*Beta(x|a,b) + e
// Used as warm-start for fit_one_context_cpp.
// Not exported to R.
// ---------------------------------------------------------------
static List fit_beta_unif_em(NumericVector x,
                               double eps_max  = 0.25,
                               int    max_iter = 150,
                               double tol      = 1e-8) {
    int n = x.size();
    NumericVector lx(n), l1mx(n);
    for (int i = 0; i < n; i++) {
        double xi = x[i];
        if (xi < 1e-8) xi = 1e-8; if (xi > 1-1e-8) xi = 1-1e-8;
        x[i] = xi; lx[i] = std::log(xi); l1mx[i] = std::log(1-xi);
    }
    double m = 0, v = 0;
    for (int i = 0; i < n; i++) m += x[i]; m /= n;
    for (int i = 0; i < n; i++) v += (x[i]-m)*(x[i]-m); v /= n;
    double mv = m*(1-m)*0.999; if (v > mv) v = mv; if (v < 1e-15) v = 1e-15;
    double phi = m*(1-m)/v - 1;
    double alpha = std::max(m*phi, 0.1), beta = std::max((1-m)*phi, 0.1);
    double eps = std::min(0.05, eps_max/2), prev_ll = -1e300;
    NumericVector ru(n), rb(n);
    for (int it = 0; it < max_iter; it++) {
        double ll = 0, sru = 0;
        for (int i = 0; i < n; i++) {
            double den = (1-eps)*dbeta_c(x[i],alpha,beta)+eps;
            if (den < 1e-300) den = 1e-300;
            ru[i] = eps/den; rb[i] = 1-ru[i]; ll += std::log(den); sru += ru[i];
        }
        eps = sru/n; if (eps > eps_max) eps = eps_max;
        NumericVector r = mle_beta_nr(lx, l1mx, rb, alpha, beta);
        alpha = r["alpha"]; beta = r["beta"];
        if (std::abs(ll-prev_ll) < tol*(1+std::abs(ll))) break;
        prev_ll = ll;
    }
    return List::create(Named("alpha")=alpha, Named("beta")=beta, Named("eps")=eps);
}

// ---------------------------------------------------------------
// Joint constrained EM for one sequence context
//
// Alternates EM updates for (pos, neg) with Euclidean projection
// onto the MLRP constraint set after each step:
//   alpha_pos >= alpha_neg  AND  beta_pos <= beta_neg
// This guarantees the posterior correction is monotone.
// ---------------------------------------------------------------
// [[Rcpp::export(.fit_one_context_cpp, rng = false)]]
List fit_one_context_cpp(NumericVector pos_data, NumericVector neg_data,
                          double eps_max_pos = 0.15,
                          double eps_max_neg = 0.25,
                          int    max_iter    = 100,
                          double tol         = 1e-7) {
    int np = pos_data.size(), nn = neg_data.size();
    NumericVector lxp(np), l1mxp(np), lxn(nn), l1mxn(nn);
    for (int i = 0; i < np; i++) {
        double xi = pos_data[i];
        if (xi<1e-8) xi=1e-8; if (xi>1-1e-8) xi=1-1e-8;
        pos_data[i]=xi; lxp[i]=std::log(xi); l1mxp[i]=std::log(1-xi);
    }
    for (int i = 0; i < nn; i++) {
        double xi = neg_data[i];
        if (xi<1e-8) xi=1e-8; if (xi>1-1e-8) xi=1-1e-8;
        neg_data[i]=xi; lxn[i]=std::log(xi); l1mxn[i]=std::log(1-xi);
    }
    List fp = fit_beta_unif_em(pos_data, eps_max_pos);
    List fn = fit_beta_unif_em(neg_data, eps_max_neg);
    double ap=fp["alpha"], bp=fp["beta"], ep=fp["eps"];
    double an=fn["alpha"], bn=fn["beta"], en=fn["eps"];
    auto proj = [](double& ap, double& bp, double& an, double& bn) {
        if (ap < an) { double m=(ap+an)/2; ap=m+1e-6; an=m-1e-6; }
        if (bp > bn) { double m=(bp+bn)/2; bp=m-1e-6; bn=m+1e-6; }
        if (ap<1e-6) ap=1e-6; if (bp<1e-6) bp=1e-6;
        if (an<1e-6) an=1e-6; if (bn<1e-6) bn=1e-6;
    };
    proj(ap, bp, an, bn);
    NumericVector rup(np), rbp(np), run(nn), rbn(nn);
    double prev_ll = -1e300; bool conv = false;
    for (int it = 0; it < max_iter; it++) {
        double ll = 0, srup = 0, srun = 0;
        for (int i = 0; i < np; i++) {
            double den = (1-ep)*dbeta_c(pos_data[i],ap,bp)+ep;
            if (den<1e-300) den=1e-300;
            rup[i]=ep/den; rbp[i]=1-rup[i]; ll+=std::log(den); srup+=rup[i];
        }
        ep = srup/np; if (ep>eps_max_pos) ep=eps_max_pos;
        NumericVector pp = mle_beta_nr(lxp, l1mxp, rbp, ap, bp);
        for (int i = 0; i < nn; i++) {
            double den = (1-en)*dbeta_c(neg_data[i],an,bn)+en;
            if (den<1e-300) den=1e-300;
            run[i]=en/den; rbn[i]=1-run[i]; ll+=std::log(den); srun+=run[i];
        }
        en = srun/nn; if (en>eps_max_neg) en=eps_max_neg;
        NumericVector pn = mle_beta_nr(lxn, l1mxn, rbn, an, bn);
        ap=pp["alpha"]; bp=pp["beta"]; an=pn["alpha"]; bn=pn["beta"];
        proj(ap, bp, an, bn);
        if (std::abs(ll-prev_ll) < tol*(1+std::abs(ll))) { conv=true; break; }
        prev_ll = ll;
    }
    return List::create(
        Named("alpha_pos")=ap, Named("beta_pos")=bp, Named("eps_pos")=ep,
        Named("alpha_neg")=an, Named("beta_neg")=bn, Named("eps_neg")=en,
        Named("converged")=conv, Named("mlrp_ok")=(ap>=an && bp<=bn));
}

// ---------------------------------------------------------------
// Pool-adjacent-violators (PAV) isotonic regression
// Modifies y[start..end) in-place to be non-decreasing.
// ---------------------------------------------------------------
static void pav_block(std::vector<double>& y, int start, int end) {
    if (end - start <= 1) return;
    struct Block { double sum; int count; };
    std::vector<Block> stk;
    stk.reserve(end - start);
    for (int i = start; i < end; i++) {
        Block b = {y[i], 1};
        while (!stk.empty() &&
               b.sum / b.count < stk.back().sum / stk.back().count) {
            b.sum   += stk.back().sum;
            b.count += stk.back().count;
            stk.pop_back();
        }
        stk.push_back(b);
    }
    int idx = start;
    for (const auto& b : stk) {
        double mean = b.sum / b.count;
        for (int j = 0; j < b.count; j++) y[idx++] = mean;
    }
}

// ---------------------------------------------------------------
// Combined posterior computation + per-context isotonic correction
//
//   1. Single loop: clamp raw_p, compute Beta-Uniform posterior.
//   2. Sort valid observations by (context, raw_p) — one std::sort.
//   3. Apply PAV per context block — O(N) total.
//   4. Scatter corrected posteriors back to output.
//
// param_idx: 1-based row indices into the params vectors.
//   NA (IntegerVector::is_na) skips the observation (output stays NA_REAL).
// ---------------------------------------------------------------
// [[Rcpp::export(.correct_posterior_cpp, rng = false)]]
NumericVector correct_posterior_cpp(
        const NumericVector& raw_p,
        const IntegerVector& param_idx,
        const NumericVector& alpha_pos,
        const NumericVector& beta_pos,
        const NumericVector& eps_pos,
        const NumericVector& alpha_neg,
        const NumericVector& beta_neg,
        const NumericVector& eps_neg,
        double pi_pos  = 0.5,
        bool   isotonic = true) {

    int N = raw_p.size();
    NumericVector out(N, NA_REAL);

    std::vector<int>    orig(N);
    std::vector<int>    ctx(N);
    std::vector<double> p_v(N);
    std::vector<double> post(N);
    int M = 0;

    for (int i = 0; i < N; i++) {
        if (IntegerVector::is_na(param_idx[i]) || !R_finite(raw_p[i])) continue;
        int c = param_idx[i] - 1;
        double p = raw_p[i];
        if (p < 1e-9)   p = 1e-9;
        if (p > 1-1e-9) p = 1-1e-9;
        double ap=alpha_pos[c], bp=beta_pos[c], ep=eps_pos[c];
        double an=alpha_neg[c], bn=beta_neg[c], en=eps_neg[c];
        double lp_ = (1.0-ep)*dbeta_c(p,ap,bp) + ep;
        double ln_ = (1.0-en)*dbeta_c(p,an,bn) + en;
        orig[M] = i; ctx[M] = c; p_v[M] = p;
        post[M] = (pi_pos*lp_) / (pi_pos*lp_ + (1.0-pi_pos)*ln_);
        M++;
    }

    if (M == 0) return out;

    std::vector<int> ord(M);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        return ctx[a] != ctx[b] ? ctx[a] < ctx[b] : p_v[a] < p_v[b];
    });

    std::vector<double> sorted_post(M);
    for (int i = 0; i < M; i++) sorted_post[i] = post[ord[i]];

    if (isotonic) {
        int start = 0;
        while (start < M) {
            int c = ctx[ord[start]];
            int end = start + 1;
            while (end < M && ctx[ord[end]] == c) end++;
            pav_block(sorted_post, start, end);
            start = end;
        }
    }

    for (int i = 0; i < M; i++) {
        double v = sorted_post[i];
        if (v < 0.0) v = 0.0; if (v > 1.0) v = 1.0;
        out[orig[ord[i]]] = v;
    }

    return out;
}
