"""Float64 Python port of the Julia design of src/policies.jl (state pass with level sums
over colex-ranked levels, count pass, Dinkelbach), used only to validate the design against
policies_reference.csv before writing the Julia code.
Usage: python3 policies_port.py <data/systems dir> <policies_reference.csv>
"""
import sys, csv, importlib.util
from math import comb
spec = importlib.util.spec_from_file_location("ref", __file__.replace("policies_port.py", "thresholds_reference.py"))
src = open(spec.origin).read().split("\ndef main():")[0]
ref = {}
exec(src, ref)


def next_comb(x):
    c = x & -x; r = x + c
    return (((r ^ x) >> 2) // c) | r


def state_pass(n, phi, q, Psi, costs, g):
    c_set, c_cmp, c_sys = costs
    fc = lambda f: c_set + c_cmp * f + c_sys
    prevC = [[fc(n)]]; prevL = [[0.0]]          # prev[col][j]
    for w in range(1, n + 1):
        f = n - w
        m = comb(n, w)
        curC = []; curL = []
        x = (1 << w) - 1
        for idx in range(m):
            if not phi[x]:
                curC.append([comb(w, j) * fc(f + j) for j in range(w + 1)])
                curL.append([0.0] * (w + 1))
            else:
                pos = [b for b in range(n) if (x >> b) & 1]
                C = [0.0] * (w + 1); L = [0.0] * (w + 1)
                for r in range(w):
                    son = [p for k, p in enumerate(pos) if k != r]
                    col = sum(comb(p, i + 1) for i, p in enumerate(son))      # colex rank (0-based)
                    for j in range(1, w + 1):
                        C[j] += prevC[col][j - 1]; L[j] += prevL[col][j - 1]
                sc = sl = 0.0
                for j in range(1, w + 1):
                    C[j] /= j; L[j] /= j
                    sc += q[w][j] / Psi[w] * C[j]; sl += q[w][j] / Psi[w] * L[j]
                cc, lc = sc, 1 / Psi[w] + sl
                call = f >= 1 and (c_set + c_cmp * f) <= cc - g * lc
                C[0] = c_set + c_cmp * f if call else cc
                L[0] = 0.0 if call else lc
                curC.append(C); curL.append(L)
            if idx < m - 1:
                x = next_comb(x)
        prevC, prevL = curC, curL
    return prevC[0][0], prevL[0][0]


def dinkel(pass_, g):
    for it in range(100):
        c, l = pass_(g)
        if c - g * l >= -1e-12 * max(1, abs(c)):
            return c / l
        g = c / l
    raise RuntimeError


systems = ref["load"](sys.argv[1])
rows = list(csv.DictReader(l for l in open(sys.argv[2]) if not l.startswith("#")))
worst = 0
cache = {}
for r in rows:
    sid = int(r["system"]); s = systems[sid]; n = len(s["edges"])
    if sid not in cache:
        cache[sid] = ref["phi_table"](s["edges"], s["directed"])
    phi = cache[sid]
    lam, Psi = ref["rates"](r["model"], r["beta"], n)
    costs = (float(r["setup"]), float(r["ratio"]) / n, 1.0)
    gs = dinkel(lambda g: state_pass(n, phi, lam, Psi, costs, g), float(r["ltmc_best_r"]))
    e = abs(gs - float(r["ltmc_state"])) / float(r["ltmc_state"])
    worst = max(worst, e)
    assert e < 1e-11, (r, gs)
print(len(rows), "cases; max relative difference of the state-based optimum:", worst)
