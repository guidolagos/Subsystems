"""Independent reference counts W_k for the 22 systems of data/systems/yeh2024_raw.txt.

W_k = number of working states with exactly k failed components. Not used by the
Julia code; it generates test/reference/yeh2024_W_reference.csv, which Pkg.test()
and scripts/check_signatures.jl compare against.

Sources of each row:
  * systems 1-5, 8-20 (undirected): the C++ reference of the Simple-delayed repo
    (test/reference/yeh2024_W_reference.csv there), which used the same edge lists;
  * systems 6D and 7D (arcs as listed in the raw file) and 22: computed here, by
    vectorized reachability over all 2^n states (numpy). This implementation is first
    checked against the C++ rows of every undirected system with n <= 21;
  * system 21: equal to system 13, since the two networks are isomorphic with the
    terminals fixed (checked here with networkx); 2^30 states are too many for numpy here.

Usage: python3 W_reference.py <raw file> <C++ reference csv of Simple-delayed>  > yeh2024_W_reference.csv
"""
import re, sys
import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism as iso


def parse_raw(path):
    out = {}
    cur = None
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        m = re.match(r'^(\d+)(D?):$', line)
        if m:
            cur = int(m[1])
            out[cur] = (m[2] == 'D', [])
            continue
        m = re.match(r'^\{"([^"]+)",\s*"([^"]+)"\},?$', line)
        out[cur][1].append((m[1], m[2]))
    return out


def counts(edges, directed, s='Start', t='End'):
    nodes = sorted({u for e in edges for u in e})
    ix = {v: i for i, v in enumerate(nodes)}
    n = len(edges)
    x = np.arange(1 << n, dtype=np.uint64)
    reach = np.full(1 << n, 1 << ix[s], dtype=np.uint64)
    while True:
        prev = reach.copy()
        for e, (u, v) in enumerate(edges):
            on = ((x >> np.uint64(e)) & np.uint64(1)).astype(bool)
            bu, bv = np.uint64(1 << ix[u]), np.uint64(1 << ix[v])
            reach = np.where(on & ((reach & bu) != 0), reach | bv, reach)
            if not directed:
                reach = np.where(on & ((reach & bv) != 0), reach | bu, reach)
        if (reach == prev).all():
            break
    works = (reach & np.uint64(1 << ix[t])) != 0
    nfailed = n - np.array([bin(i).count('1') for i in range(1 << n)])
    return [int(works[nfailed == k].sum()) for k in range(n + 1)]


def read_cpp(path):
    ref = {}
    for line in open(path):
        if line.startswith('#') or line.startswith('system') or not line.strip():
            continue
        i, k, w = map(int, line.split(','))
        ref.setdefault(i, []).append(w)
    return ref


def graph(edges):
    g = nx.MultiGraph()
    g.add_edges_from(edges)
    for v in g:
        g.nodes[v]['r'] = 's' if v == 'Start' else ('t' if v == 'End' else '')
    return g


raw = parse_raw(sys.argv[1])
cpp = read_cpp(sys.argv[2])

# 1. validate this implementation against the C++ counts (undirected systems, n <= 21)
for i in sorted(cpp):
    directed, edges = raw[i]
    if len(edges) <= 21:
        assert counts(edges, False) == cpp[i], f"python != C++ for system {i}"
        print(f"# python == C++ for system {i} (undirected)", file=sys.stderr)

# 2. system 21 is system 13 redrawn
assert nx.is_isomorphic(graph(raw[21][1]), graph(raw[13][1]),
                        node_match=iso.categorical_node_match('r', ''))

W = {}
for i, (directed, edges) in sorted(raw.items()):
    if i in (6, 7, 22):
        W[i] = counts(edges, directed)
    elif i == 21:
        W[i] = cpp[13]
    else:
        assert not directed
        W[i] = cpp[i]
    assert len(W[i]) == len(edges) + 1 and W[i][0] == 1 and W[i][-1] == 0

print("# W[k] = number of working states with exactly k failed components, for the systems of")
print("# data/systems/yeh2024_*.csv, computed independently of the Julia code by W_reference.py:")
print("# systems 1-5 and 8-20 from the C++ reference of the Simple-delayed repository (same edge lists,")
print("# undirected); systems 6 and 7 (directed, arcs as in the raw file) and 22 by numpy enumeration,")
print("# validated against the C++ rows; system 21 = system 13 (isomorphic networks).")
print("system,k,W")
for i in sorted(W):
    for k, w in enumerate(W[i]):
        print(f"{i},{k},{w}")
