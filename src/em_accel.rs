//! Fixed-point acceleration used by the abundance EM.
//!
//! Both accelerators treat the existing M-step as an opaque map `F(x)`.  The
//! expensive map may be parallel; all extrapolation reductions are deliberately
//! sequential so their result does not depend on rayon scheduling.

const DAAREM_ORDER: usize = 5;
const DAAREM_A1: f64 = 1.2;
const DAAREM_KAPPA: f64 = 25.0;
const DAAREM_RHO: f64 = 0.95;

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn feasible(v: &[f64]) -> bool {
    v.iter().all(|x| x.is_finite() && *x >= 0.0)
}

/// SQUAREM SqS3. `step` returns the convergence distance for a completed map.
#[allow(clippy::too_many_arguments)]
pub(crate) fn squarem(
    f: &mut impl FnMut(&[f64], &mut [f64]),
    x: &mut Vec<f64>,
    scratch: &mut Vec<f64>,
    max_eval: u32,
    min_eval: u32,
    evals: &mut u32,
    tolerance: f64,
    distance: impl Fn(&[f64], &[f64]) -> f64,
) -> bool {
    let n = x.len();
    let mut x1 = vec![0.0; n];
    let mut x2 = vec![0.0; n];
    let mut r = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut proposal = vec![0.0; n];

    while *evals < max_eval {
        // A partial cycle at the end of the budget is ordinary EM.
        if max_eval - *evals < 3 {
            f(x, scratch);
            *evals += 1;
            let d = distance(x, scratch);
            std::mem::swap(x, scratch);
            if *evals >= min_eval && d.is_finite() && d < tolerance {
                return true;
            }
            continue;
        }

        f(x, &mut x1);
        f(&x1, &mut x2);
        *evals += 2;

        let mut r2 = 0.0;
        let mut v2 = 0.0;
        for i in 0..n {
            r[i] = x1[i] - x[i];
            v[i] = (x2[i] - x1[i]) - r[i];
            r2 += r[i] * r[i];
            v2 += v[i] * v[i];
        }
        if !r2.is_finite() || !v2.is_finite() || v2 <= f64::EPSILON {
            let d = distance(x, &x2);
            std::mem::swap(x, &mut x2);
            if *evals >= min_eval && d.is_finite() && d < tolerance {
                return true;
            }
            continue;
        }

        let mut alpha = -(r2 / v2).sqrt().min(1.0e6);
        alpha = alpha.min(-1.0);
        let mut accepted = false;
        for _ in 0..=30 {
            let a2 = alpha * alpha;
            for i in 0..n {
                proposal[i] = x[i] - 2.0 * alpha * r[i] + a2 * v[i];
            }
            if feasible(&proposal) {
                accepted = true;
                break;
            }
            alpha = (alpha - 1.0) / 2.0;
        }
        if !accepted {
            proposal.copy_from_slice(&x2);
        }

        f(&proposal, scratch);
        *evals += 1;
        let d = distance(x, scratch);
        std::mem::swap(x, scratch);
        if *evals >= min_eval && d.is_finite() && d < tolerance {
            return true;
        }
    }
    false
}

fn jacobi_eig(a: &mut [f64], n: usize, eval: &mut [f64], evec: &mut [f64]) {
    for r in 0..n {
        for c in 0..n {
            evec[r * n + c] = f64::from(r == c);
        }
    }
    for _ in 0..50 {
        let off = (0..n)
            .flat_map(|p| ((p + 1)..n).map(move |q| (p, q)))
            .map(|(p, q)| a[p * n + q].powi(2))
            .sum::<f64>();
        if off <= 1e-30 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                if apq.abs() <= 1e-300 {
                    continue;
                }
                let phi = 0.5 * (2.0 * apq).atan2(a[q * n + q] - a[p * n + p]);
                let (s, c) = phi.sin_cos();
                for k in 0..n {
                    let (akp, akq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (evec[k * n + p], evec[k * n + q]);
                    evec[k * n + p] = c * vkp - s * vkq;
                    evec[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    for i in 0..n {
        eval[i] = a[i * n + i];
    }
}

#[allow(clippy::too_many_arguments)]
fn damping_find(
    uy_sq: &[f64],
    dvec: &[f64],
    sk: f64,
    ftf: f64,
    lambda_start: f64,
    r_start: f64,
) -> (f64, f64) {
    let valid: Vec<usize> = (0..dvec.len()).filter(|&i| dvec[i] > 0.0).collect();
    if valid.is_empty() {
        return (lambda_start, r_start);
    }
    let pow = DAAREM_KAPPA - sk;
    let target = (-0.5 * (1.0 + DAAREM_A1.powf(pow)).ln()).exp();
    let ls_norm = valid
        .iter()
        .map(|&i| uy_sq[i] / dvec[i].powi(2))
        .sum::<f64>()
        .sqrt();
    let vk = target * ls_norm;
    if vk == 0.0 {
        return (lambda_start, r_start);
    }
    let mut lambda = lambda_start - r_start / vk;
    let denom = valid
        .iter()
        .map(|&i| uy_sq[i] / dvec[i].powi(4))
        .sum::<f64>();
    let mut lower = ls_norm * (ls_norm - vk) / denom;
    let mut upper = ftf / vk;
    let low_stop = (-0.5 * (1.0 + DAAREM_A1.powf(pow + 0.5)).ln()).exp();
    let high_stop = (-0.5 * (1.0 + DAAREM_A1.powf(pow - 0.5)).ln()).exp();
    let mut ratio = 0.0;
    let mut s_norm = 0.0;
    for _ in 0..10 {
        if lambda <= lower || lambda >= upper {
            lambda = (0.0001 * upper).max((lower * upper).sqrt());
        }
        let mut sn = 0.0;
        let mut der = 0.0;
        for &i in &valid {
            let denom = dvec[i].powi(2) + lambda;
            let dl = (dvec[i] / denom).powi(2);
            sn += uy_sq[i] * dl;
            der += uy_sq[i] * dl / denom;
        }
        s_norm = sn.sqrt();
        let phi = s_norm - vk;
        ratio = phi / (-der / s_norm);
        if s_norm <= high_stop * ls_norm && s_norm >= low_stop * ls_norm {
            break;
        }
        if phi < 0.0 {
            upper = lambda;
        }
        lower = lower.max(lambda - ratio);
        lambda -= s_norm * ratio / vk;
    }
    (lambda, s_norm * ratio)
}

/// Objective-free damped Anderson acceleration with history restarts and
/// residual-monotonicity control (DAAREM).
#[allow(clippy::too_many_arguments, clippy::ptr_arg)]
pub(crate) fn daarem(
    f: &mut impl FnMut(&[f64], &mut [f64]),
    x: &mut Vec<f64>,
    scratch: &mut Vec<f64>,
    max_eval: u32,
    min_eval: u32,
    evals: &mut u32,
    tolerance: f64,
    distance: impl Fn(&[f64], &[f64]) -> f64,
) -> bool {
    let n = x.len();
    if n == 0 {
        return true;
    }
    let order = DAAREM_ORDER.min(n.div_ceil(2)).max(1);
    if max_eval - *evals < 2 {
        return false;
    }
    let mut xold = x.clone();
    let mut xnew = vec![0.0; n];
    f(&xold, &mut xnew);
    *evals += 1;
    let mut fold: Vec<f64> = (0..n).map(|i| xnew[i] - xold[i]).collect();
    f(&xnew, scratch);
    *evals += 1;
    let mut fnew: Vec<f64> = (0..n).map(|i| scratch[i] - xnew[i]).collect();
    let mut residual_norm = norm(&fnew);
    let mut fdiff = vec![vec![0.0; n]; order];
    let mut xdiff = vec![vec![0.0; n]; order];
    let mut proposal = vec![0.0; n];
    let mut mapped = vec![0.0; n];
    let mut count = 0usize;
    let mut accepted = 0.0;
    let mut lambda = 100_000.0;
    let mut penalty = 0.0;
    let mut cycle = 1i32;

    while *evals < max_eval {
        for i in 0..n {
            fdiff[count][i] = fnew[i] - fold[i];
            xdiff[count][i] = xnew[i] - xold[i];
        }
        count += 1;
        let m = count;
        let mut gram = vec![0.0; m * m];
        let mut b = vec![0.0; m];
        for a in 0..m {
            for c in a..m {
                let value = (0..n).map(|i| fdiff[a][i] * fdiff[c][i]).sum();
                gram[a * m + c] = value;
                gram[c * m + a] = value;
            }
            b[a] = (0..n).map(|i| fdiff[a][i] * fnew[i]).sum();
        }
        let mut eigenvalues = vec![0.0; m];
        let mut eigenvectors = vec![0.0; m * m];
        jacobi_eig(&mut gram, m, &mut eigenvalues, &mut eigenvectors);
        let singular: Vec<f64> = eigenvalues.iter().map(|v| v.max(0.0).sqrt()).collect();
        let mut uy = vec![0.0; m];
        let mut ftf_sq = 0.0;
        for a in 0..m {
            let vtb = (0..m).map(|r| eigenvectors[r * m + a] * b[r]).sum::<f64>();
            ftf_sq += vtb * vtb;
            uy[a] = if singular[a] > 0.0 {
                vtb / singular[a]
            } else {
                0.0
            };
        }
        let uy_sq: Vec<f64> = uy.iter().map(|v| v * v).collect();
        (lambda, penalty) =
            damping_find(&uy_sq, &singular, accepted, ftf_sq.sqrt(), lambda, penalty);
        let dd: Vec<f64> = (0..m)
            .map(|a| singular[a] * uy[a] / (singular[a].powi(2) + lambda))
            .collect();
        let gamma: Vec<f64> = (0..m)
            .map(|r| (0..m).map(|a| eigenvectors[r * m + a] * dd[a]).sum())
            .collect();
        for i in 0..n {
            let xg = (0..m).map(|a| xdiff[a][i] * gamma[a]).sum::<f64>();
            let fg = (0..m).map(|a| fdiff[a][i] * gamma[a]).sum::<f64>();
            proposal[i] = (xnew[i] - xg) + (fnew[i] - fg);
        }

        let old = xnew.clone();
        let use_proposal = feasible(&proposal);
        if use_proposal {
            f(&proposal, &mut mapped);
            *evals += 1;
        }
        let proposed_residual: Vec<f64> = (0..n).map(|i| mapped[i] - proposal[i]).collect();
        let proposed_norm = norm(&proposed_residual);
        std::mem::swap(&mut xold, &mut xnew);
        if use_proposal
            && proposed_norm.is_finite()
            && proposed_norm <= residual_norm * (1.0 + DAAREM_RHO.powi(cycle))
        {
            std::mem::swap(&mut fold, &mut fnew);
            xnew.copy_from_slice(&proposal);
            fnew.copy_from_slice(&proposed_residual);
            accepted += 1.0;
            residual_norm = proposed_norm;
        } else {
            std::mem::swap(&mut fold, &mut fnew);
            for i in 0..n {
                xnew[i] = xold[i] + fold[i];
            }
            if *evals < max_eval {
                f(&xnew, &mut mapped);
                *evals += 1;
                for i in 0..n {
                    fnew[i] = mapped[i] - xnew[i];
                }
                residual_norm = norm(&fnew);
            }
        }
        let d = distance(&old, &xnew);
        if *evals >= min_eval && d.is_finite() && d < tolerance {
            x.copy_from_slice(&xnew);
            return true;
        }
        if count == order {
            count = 0;
        }
        cycle += 1;
    }
    x.copy_from_slice(&xnew);
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    fn distance(a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max)
    }

    fn slow_map(src: &[f64], dst: &mut [f64]) {
        dst[0] = 0.98 * src[0] + 0.02;
        dst[1] = 2.0 - dst[0];
    }

    #[test]
    fn squarem_reaches_nonnegative_fixed_point() {
        let mut x = vec![0.0, 2.0];
        let mut scratch = vec![0.0; 2];
        let mut evals = 0;
        let converged = squarem(
            &mut slow_map,
            &mut x,
            &mut scratch,
            500,
            1,
            &mut evals,
            1e-8,
            distance,
        );
        assert!(converged);
        assert!((x[0] - 1.0).abs() < 1e-6, "{x:?}");
        assert!(x.iter().all(|v| *v >= 0.0));
        assert!(evals <= 500);
    }

    #[test]
    fn daarem_reaches_nonnegative_fixed_point() {
        let mut x = vec![0.0, 2.0];
        let mut scratch = vec![0.0; 2];
        let mut evals = 0;
        let converged = daarem(
            &mut slow_map,
            &mut x,
            &mut scratch,
            500,
            1,
            &mut evals,
            1e-8,
            distance,
        );
        assert!(converged);
        assert!((x[0] - 1.0).abs() < 1e-6, "{x:?}");
        assert!(x.iter().all(|v| *v >= 0.0));
        assert!(evals <= 500);
    }
}
