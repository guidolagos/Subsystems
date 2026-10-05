"""Independent reference for the MTTF-threshold policies pi^gamma of the paper.

Follows the paper literally, without signatures:
  * Psi of each model, calibrated to Psi(1) = 1 and 2 - Psi(2) = beta, in mpmath (60 digits);
  * rate from a state with w working components to ONE given sub-state with j fewer working
    components: lambda^(w)_j = sum_{l=w-j}^{w} (-1)^(l-(w-j)+1) binom(j, w-l) Psi(l)
    (eq. def:lambda of the companion paper), in mpmath;
  * MTTF at every state by the recursion eq:mttf over the 2^n states;
  * for a threshold gamma, the policy repairs at the states with m(x) < gamma m(1); its
    mean time to repair T, mean number of failed components at repair N and probability that
    the repair is corrective (system failure) P are obtained from the linear systems
    sys:m rep and sys:c rep, solved by recursion over the states (sub-states first).
The breakpoints are the distinct values m(x)/m(1) > 0, in decreasing order (ties within
1e-12 m(1) merged); breakpoint i means: continue at the states with m(x) >= gamma_i m(1).

Usage: python3 thresholds_reference.py <data/systems dir> > thresholds_reference.csv
"""
import sys, csv, math
import numpy as np
import mpmath as mp

mp.mp.dps = 60


def load(dirpath):
    systems = {}
    for r in csv.DictReader(open(f"{dirpath}/yeh2024_systems.csv")):
        systems[int(r["system"])] = dict(directed=r["directed"] == "true", edges=[])
    for r in csv.DictReader(open(f"{dirpath}/yeh2024_edges.csv")):
        systems[int(r["system"])]["edges"].append((r["from"], r["to"]))
    return systems


def phi_table(edges, directed):
    n = len(edges)
    phi = np.zeros(1 << n, dtype=bool)
    for x in range(1 << n):
        adj = {}
        for e, (u, v) in enumerate(edges):
            if (x >> e) & 1:
                adj.setdefault(u, []).append(v)
                if not directed:
                    adj.setdefault(v, []).append(u)
        seen, stack = {"Start"}, ["Start"]
        while stack:
            u = stack.pop()
            for v in adj.get(u, []):
                if v not in seen:
                    seen.add(v); stack.append(v)
        phi[x] = "End" in seen
    return phi


def psi_fun(name, a):
    if name == "Exp":     return lambda x: mp.mpf(x)
    if name == "LFMO-C":  return lambda x: x / (1 + a * (x - 1))
    if name == "LFMO-G":  return lambda x: mp.log(1 + x / a) / mp.log(1 + 1 / a)
    if name == "LFMO-IG": return lambda x: (mp.sqrt(2 * x + a**2) - a) / (mp.sqrt(2 + a**2) - a)
    if name == "LFMO-S":  return lambda x: mp.mpf(x)**a
    if name == "BFR":     return lambda x: (1 - (1 - a)**x) / a
    if name == "LSh":     return lambda x: (1 - a) * x + a
    raise ValueError(name)


def calibrate(name, beta):
    b = mp.mpf(beta)
    if name == "Exp": return mp.mpf(0)
    if name == "LFMO-C": return b / (2 - b)
    if name in ("BFR", "LSh"): return b
    if name == "LFMO-S": return mp.log(2 - b, 2)
    # LFMO-G, LFMO-IG: beta(eta) = 2 - Psi(2) decreases in eta; bisection on log(eta)
    g = lambda le: 2 - psi_fun(name, mp.e**le)(2) - b
    lo, hi = mp.mpf(-60), mp.mpf(60)
    assert g(lo) > 0 > g(hi)
    for _ in range(400):
        mid = (lo + hi) / 2
        if g(mid) > 0: lo = mid
        else: hi = mid
    return mp.e**((lo + hi) / 2)


def rates(name, beta, n):
    a = calibrate(name, beta)
    P = psi_fun(name, a)
    Psi = [mp.mpf(0)] + [P(x) for x in range(1, n + 1)]
    assert abs(Psi[1] - 1) < mp.mpf(10)**-40 and abs(2 - Psi[2] - mp.mpf(beta)) < mp.mpf(10)**-40
    lam = np.zeros((n + 1, n + 1))   # lam[w, j]
    for w in range(1, n + 1):
        tot = mp.mpf(0)
        for j in range(1, w + 1):
            s = mp.mpf(0)
            for l in range(w - j, w + 1):
                s += (-1)**(l - (w - j) + 1) * mp.binomial(j, w - l) * Psi[l]
            assert s > -mp.mpf(10)**-40
            s = max(s, mp.mpf(0))
            lam[w, j] = float(s)
            tot += mp.binomial(w, j) * s
        assert abs(tot - Psi[w]) < mp.mpf(10)**-40 * Psi[w]
    return lam, [float(p) for p in Psi]


def subsets(x):
    """Proper sub-states y < x."""
    y = (x - 1) & x
    while True:
        yield y
        if y == 0: break
        y = (y - 1) & x


def run(n, phi, lam, Psi, idx_select):
    order = sorted(range(1 << n), key=lambda x: bin(x).count("1"))
    pc = [bin(x).count("1") for x in range(1 << n)]
    m = np.zeros(1 << n)
    for x in order:
        if not phi[x]: continue
        w = pc[x]
        s = 1.0
        for y in subsets(x):
            s += lam[w, w - pc[y]] * m[y]
        m[x] = s / Psi[w]
    m1 = m[(1 << n) - 1]
    vals = sorted({v for v in m if v > 0}, reverse=True)
    groups = []
    for v in vals:
        if groups and groups[-1] - v <= 1e-12 * m1: continue
        groups.append(v)
    G = len(groups)
    sel = sorted(set(i for i in idx_select(G)))
    thr = np.array([groups[i] for i in sel]) - 1e-12 * m1
    k = len(sel)
    T = np.zeros((1 << n, k)); N = np.zeros((1 << n, k)); Pf = np.zeros((1 << n, k))
    for x in order:
        if not phi[x]: continue
        w = pc[x]
        cont = m[x] >= thr                       # bool vector over the selected gammas
        t = np.full(k, 1.0); nn = np.zeros(k); pp = np.zeros(k)
        for y in subsets(x):
            r = lam[w, w - pc[y]]
            if r == 0: continue
            if phi[y]:
                cy = m[y] >= thr
                t += r * np.where(cy, T[y], 0.0)
                nn += r * np.where(cy, N[y], n - pc[y])
                pp += r * np.where(cy, Pf[y], 0.0)
            else:
                nn += r * (n - pc[y])
                pp += r
        T[x] = np.where(cont, t / Psi[w], 0.0)
        N[x] = np.where(cont, nn / Psi[w], n - w)
        Pf[x] = np.where(cont, pp / Psi[w], 0.0)
    full = (1 << n) - 1
    mps = sum(1 for x in range(1 << n) if phi[x] and all(not phi[x ^ (1 << e)] for e in range(n) if (x >> e) & 1))
    return m1, G, sel, [groups[i] / m1 for i in sel], T[full], N[full], Pf[full], mps


def main():
    systems = load(sys.argv[1])
    models = [("Exp", "0"), ("LFMO-C", "0.2"), ("LFMO-G", "0.05"), ("LFMO-IG", "0.5"),
              ("LFMO-S", "0.1"), ("BFR", "0.5"), ("LSh", "0.2")]
    cases = [(1, models), (3, models), (9, models), (6, [models[0], models[1], models[6]])]
    out = csv.writer(sys.stdout, lineterminator="\n")
    print("# Generated by test/reference/thresholds_reference.py (independent of the Julia code).")
    out.writerow(["system", "model", "beta", "m1", "n_breakpoints", "mps", "index", "gamma", "T", "N", "P"])
    for sid, mods in cases:
        s = systems[sid]
        n = len(s["edges"])
        phi = phi_table(s["edges"], s["directed"])
        for name, beta in mods:
            lam, Psi = rates(name, beta, n)
            if n <= 8:
                sel = lambda G: range(G)
            else:
                sel = lambda G: [0, 1, 2, 3, G // 4, G // 2, (3 * G) // 4, G - 3, G - 2, G - 1]
            m1, G, idx, gam, T, N, P, mps = run(n, phi, lam, Psi, sel)
            # identities: gamma -> 0+ is corrective; T = m(1), P = 1
            assert abs(T[-1] - m1) < 1e-12 * m1 and abs(P[-1] - 1) < 1e-12
            for i, g, t, nn, p in zip(idx, gam, T, N, P):
                out.writerow([sid, name, beta, repr(float(m1)), G, mps, i + 1, repr(float(g)), repr(float(t)), repr(float(nn)), repr(float(p))])
            print(f"# system {sid} {name} {beta}: n = {n}, {G} breakpoints, m1 = {m1:.12g}", file=sys.stderr)


main()
