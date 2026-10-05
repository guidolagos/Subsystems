"""Independent reference for the benchmark policies of src/policies.jl (no delay).

Follows the paper literally ("Optimal policy" in Section sec:policies), with explicit
sums over all sub-states x' < x (no level sums, no signatures):
  * a policy is evaluated by solving the linear systems sys:m rep and sys:c rep by
    recursion over the states (sub-states first); repair cost c_set + f c_cmp, plus
    c_sys at failed states (f = number of failed components);
  * the optimal state-based policy by the algorithm of the paper (Dinkelbach / policy
    iteration: h-recursion of step 2, evaluation of step 3), and, for the bridge
    (system 1), also by brute force over all 2^15 state-based policies;
  * the optimal count-based policy by brute force over all sets S of counts f >= 1 at
    which the policy repairs, each evaluated with the survivor-count chain; for
    systems 1 and 3 every count policy is also evaluated with the explicit state
    equations (checks the count-chain evaluation);
  * r-out-of-n policies (repair when f >= r or the system failed), r = 1..n, with the
    explicit state equations.
Rates from thresholds_reference.py (Psi calibrated in 60-digit arithmetic).

Usage: python3 policies_reference.py <data/systems dir> > policies_reference.csv
"""
import sys, csv, itertools, importlib.util
import numpy as np

spec = importlib.util.spec_from_file_location("ref", __file__.replace("policies_reference.py", "thresholds_reference.py"))
src = open(spec.origin).read().split("\ndef main():")[0]
ref = {}
exec(src, ref)          # load, phi_table, rates (Psi, lambda in mpmath)


class States:
    def __init__(self, edges, directed):
        self.n = n = len(edges)
        self.N = 1 << n
        self.phi = ref["phi_table"](edges, directed)
        self.pc = np.array([bin(x).count("1") for x in range(self.N)])
        self.order = sorted(range(self.N), key=lambda x: self.pc[x])
        self.subs = {}
        for x in range(self.N):
            if not self.phi[x]:
                continue
            bits = [b for b in range(n) if (x >> b) & 1]
            idx = np.arange(1 << len(bits), dtype=np.int64)
            y = np.zeros(1 << len(bits), dtype=np.int64)
            for k, b in enumerate(bits):
                y |= ((idx >> k) & 1) << b
            y = y[y != x]
            self.subs[x] = (y, self.pc[x] - self.pc[y])        # sub-states and j = |x| - |y|


def evaluate(S, lam, Psi, costs, repair):
    """LTMC of the policy that repairs at working x iff repair(x) (never at the full state)."""
    n = S.n
    c_set, c_cmp, c_sys = costs
    m = np.zeros(S.N); c = np.zeros(S.N)
    full = S.N - 1
    for x in S.order:
        f = n - S.pc[x]
        if not S.phi[x]:
            c[x] = c_set + c_cmp * f + c_sys
            continue
        if x != full and repair(x):
            c[x] = c_set + c_cmp * f
            continue
        w = S.pc[x]
        y, j = S.subs[x]
        p = lam[w, j] / Psi[w]
        m[x] = 1 / Psi[w] + np.dot(p, m[y])
        c[x] = np.dot(p, c[y])
    return c[full] / m[full]


def optimal_state(S, lam, Psi, costs, g0):
    """The algorithm of the paper: h-recursion (step 2), evaluation by sys:m rep, sys:c rep (step 3)."""
    n = S.n
    c_set, c_cmp, c_sys = costs
    full = S.N - 1
    g = g0
    for it in range(1, 200):
        h = np.zeros(S.N)
        rep = np.zeros(S.N, dtype=bool)
        for x in S.order:
            f = n - S.pc[x]
            if not S.phi[x]:
                h[x] = c_set + c_cmp * f + c_sys
                continue
            w = S.pc[x]
            y, j = S.subs[x]
            cont = -g / Psi[w] + np.dot(lam[w, j] / Psi[w], h[y])
            if x == full:
                h[x] = cont
            else:
                K = c_set + c_cmp * f
                rep[x] = K <= cont
                h[x] = min(K, cont)
        gnew = evaluate(S, lam, Psi, costs, lambda x: rep[x])
        if h[full] >= -1e-12 * max(1.0, abs(g)):
            return gnew, rep, it
        g = gnew
    raise RuntimeError("no convergence")


def count_eval(Sbar, lam, Psi, n, costs, Sset):
    """LTMC of the count policy repairing at working states with f in Sset, via the survivor-count chain."""
    from math import comb
    c_set, c_cmp, c_sys = costs
    C = np.zeros(n + 1); T = np.zeros(n + 1)
    for f in range(n - 1, -1, -1):
        if Sbar[f] == 0:
            continue
        if f >= 1 and f in Sset:
            C[f] = c_set + c_cmp * f; T[f] = 0.0
            continue
        w = n - f
        cc = 0.0; tt = 1 / Psi[w]
        for j in range(1, w + 1):
            Pj = comb(w, j) * lam[w, j] / Psi[w]
            ps = Sbar[f + j] / Sbar[f]
            cc += Pj * (ps * C[f + j] + (1 - ps) * (c_set + c_cmp * (f + j) + c_sys))
            tt += Pj * ps * T[f + j]
        C[f] = cc; T[f] = tt
    return C[0] / T[0]


def main():
    from math import comb
    systems = ref["load"](sys.argv[1])
    allmodels = [("Exp", "0"), ("LFMO-C", "0.2"), ("LFMO-G", "0.05"), ("LFMO-IG", "0.5"),
                 ("LFMO-S", "0.1"), ("BFR", "0.5"), ("LSh", "0.2")]
    costs4 = [(0.1, 0.0), (2 / 3, 0.05), (1.0, 0.1), (10.0, 0.5)]
    cases = [(1, allmodels, costs4), (3, allmodels, costs4),
             (9, [allmodels[0], allmodels[3], allmodels[6]], costs4),
             (6, [allmodels[0], allmodels[5]], [(0.1, 0.0), (1.0, 0.1)])]
    out = csv.writer(sys.stdout, lineterminator="\n")
    print("# Generated by test/reference/policies_reference.py (independent of the Julia code).")
    out.writerow(["system", "model", "beta", "ratio", "setup", "ltmc_state", "ltmc_state_brute", "ltmc_count",
                  "best_r", "ltmc_best_r", "ltmc_r1", "ltmc_rn"])
    for sid, mods, cps in cases:
        s = systems[sid]
        S = States(s["edges"], s["directed"])
        n = S.n
        W = [int(sum(1 for x in range(S.N) if S.phi[x] and S.pc[x] == n - k)) for k in range(n + 1)]
        Sbar = [W[k] / comb(n, k) for k in range(n + 1)]
        reach = [f for f in range(1, n) if Sbar[f] > 0]
        for name, beta in mods:
            lam, Psi = ref["rates"](name, beta, n)
            for ratio, setup in cps:
                costs = (setup, ratio / n, 1.0)
                Lr = [evaluate(S, lam, Psi, costs, lambda x, r=r: n - S.pc[x] >= r) for r in range(1, n + 1)]
                ib = int(np.argmin(Lr))
                best_count = np.inf
                for k in range(len(reach) + 1):
                    for Sset in itertools.combinations(reach, k):
                        v = count_eval(Sbar, lam, Psi, n, costs, set(Sset))
                        if sid in (1, 3):
                            v2 = evaluate(S, lam, Psi, costs, lambda x, Ss=set(Sset): n - S.pc[x] in Ss)
                            assert abs(v - v2) <= 1e-12 * v2, (sid, name, Sset, v, v2)
                        best_count = min(best_count, v)
                assert best_count <= min(Lr) * (1 + 1e-12)
                gs, rep, it = optimal_state(S, lam, Psi, costs, min(Lr))
                assert gs <= best_count * (1 + 1e-12)
                brute = ""
                if sid == 1 and name in ("Exp", "LFMO-IG", "LSh"):
                    free = [x for x in range(S.N) if S.phi[x] and x != S.N - 1]
                    best = np.inf
                    for mask in range(1 << len(free)):
                        Rset = {free[i] for i in range(len(free)) if (mask >> i) & 1}
                        best = min(best, evaluate(S, lam, Psi, costs, lambda x: x in Rset))
                    assert abs(best - gs) <= 1e-12 * gs, (name, best, gs)
                    brute = repr(float(best))
                out.writerow([sid, name, beta, repr(ratio), repr(setup), repr(float(gs)), brute, repr(float(best_count)),
                              ib + 1, repr(float(Lr[ib])), repr(float(Lr[0])), repr(float(Lr[-1]))])
                print(f"# system {sid} {name} {beta} ratio {ratio:.3g} setup {setup}: state {gs:.10g} "
                      f"count {best_count:.10g} best r {ib + 1} ({Lr[ib]:.10g}); {it} iterations", file=sys.stderr)
                sys.stdout.flush()


main()
