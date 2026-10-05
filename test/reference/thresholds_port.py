"""Float64 Python port of the Julia design (src/sweep.jl, src/thresholds.jl), used only to
validate the design against thresholds_reference.py before writing the Julia code.

Design: one sweep computes the counting vector W^x (Lemma lemma:Wk) of every state and groups
the working states into classes of equal (|x|, W^x), with multiplicities. For an LFMO model:
  m(x) = sum_{k<w} W^x_k / binom(w,k) * delta^w_k,   delta^w_k = p^w_k / Psi(w-k),
where p^w_k is the probability that the survivor count of a w-component system visits w-k;
and, since m is monotone (pathwise coupling), the continuation region {m >= gamma m(1)} is an
up-set, so a state x in it is visited as in the free chain: expected time tau_w =
p^n_{n-w} / (binom(n,w) Psi(w)). Per state of the region:
  time  tau_w,
  failed components accumulated  tau_w * w * Psi(1),
  probability of a system failure at its exit  tau_w * sum_j q_w(j) (binom(w,j) - W^x_j),
with q_w(j) the rate to one given sub-state with j fewer working components.
Usage: python3 thresholds_port.py <data/systems dir> <thresholds_reference.csv>
"""
import sys, csv
from math import comb
import importlib.util
spec = importlib.util.spec_from_file_location("ref", __file__.replace("thresholds_port.py", "thresholds_reference.py"))
src = open(spec.origin).read().split("\ndef main():")[0]
ref = {}
exec(src, ref)          # reuse only: load, phi_table, rates (Psi and lambda in mpmath)


def sweep(n, phi):
    W = {}
    classes = {}
    mps = 0
    for w in range(n + 1):
        for x in (x for x in range(1 << n) if bin(x).count("1") == w):
            if not phi[x]:
                W[x] = (0,) * (w + 1); continue
            acc = [0] * w
            y = x
            while y:
                b = y & -y; y ^= b
                c = W[x ^ b]
                for j in range(w): acc[j] += c[j]
            assert all(acc[j] % (j + 1) == 0 for j in range(w))
            W[x] = (1,) + tuple(acc[j] // (j + 1) for j in range(w))
            if w >= 1 and W[x][1] == 0: mps += 1
            classes[(w, W[x])] = classes.get((w, W[x]), 0) + 1
    return classes, mps


def tables(n, lam, Psi):
    P = [[0.0] * (n + 1) for _ in range(n + 1)]
    for w in range(1, n + 1):
        for j in range(1, w + 1):
            P[w][j] = lam[w, j] * comb(w, j) / Psi[w]
    delta = [[0.0] * (n + 1) for _ in range(n + 1)]
    for w in range(1, n + 1):
        p = [0.0] * w; p[0] = 1.0
        for k in range(1, w):
            p[k] = sum(p[i] * P[w - i][k - i] for i in range(k))
        for k in range(w):
            delta[w][k] = p[k] / Psi[w - k]
    tau = [0.0] + [delta[n][n - w] / comb(n, w) for w in range(1, n + 1)]
    return delta, tau


def curve(n, classes, lam, Psi):
    delta, tau = tables(n, lam, Psi)
    rows = []
    for (w, Wx), c in classes.items():
        m = sum(Wx[k] / comb(w, k) * delta[w][k] for k in range(w))
        t = c * tau[w]
        nf = c * tau[w] * w * Psi[1]
        pf = c * tau[w] * sum(lam[w, j] * (comb(w, j) - Wx[j]) for j in range(1, w + 1))
        rows.append((m, t, nf, pf))
    rows.sort(key=lambda r: -r[0])
    m1 = rows[0][0]
    out = []           # (gamma, T, N, P) cumulative per group
    T = N = Pp = 0.0
    rep = None
    for m, t, nf, pf in rows:
        if rep is not None and rep - m > 1e-12 * m1:
            out.append((rep / m1, T, N, Pp))
        if rep is None or rep - m > 1e-12 * m1: rep = m
        T += t; N += nf; Pp += pf
    out.append((rep / m1, T, N, Pp))
    return m1, out


systems = ref["load"](sys.argv[1])
refrows = list(csv.DictReader(l for l in open(sys.argv[2]) if not l.startswith("#")))
keys = sorted({(int(r["system"]), r["model"], r["beta"]) for r in refrows}, key=lambda k: (k[0], k[1]))
worst = 0.0
cache = {}
for sid, name, beta in keys:
    s = systems[sid]; n = len(s["edges"])
    if sid not in cache:
        cache[sid] = sweep(n, ref["phi_table"](s["edges"], s["directed"]))
    classes, mps = cache[sid]
    lam, Psi = ref["rates"](name, beta, n)
    m1, out = curve(n, classes, lam, Psi)
    rr = [r for r in refrows if (int(r["system"]), r["model"], r["beta"]) == (sid, name, beta)]
    assert int(rr[0]["n_breakpoints"]) == len(out), (sid, name, len(out), rr[0]["n_breakpoints"])
    assert int(rr[0]["mps"]) == mps
    assert abs(float(rr[0]["m1"]) - m1) <= 1e-12 * m1
    for r in rr:
        g, T, N, P = out[int(r["index"]) - 1]
        for a, b in ((g, r["gamma"]), (T, r["T"]), (N, r["N"]), (P, r["P"])):
            b = float(b); err = abs(a - b) / max(abs(b), 1e-300) if b != 0 else abs(a)
            worst = max(worst, err)
            assert err < 1e-10, (sid, name, beta, r["index"], a, b)
    print(f"system {sid} {name} {beta}: {len(classes)} classes, {len(out)} breakpoints ok")
print("max relative difference:", worst)
