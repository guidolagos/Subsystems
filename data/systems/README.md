# Test systems: Yeh et al. (2024)

Source: W.-C. Yeh, G. Yang, C.-L. Huang, "A New Node-Based Algorithm for
Identifying the Complete Minimal Cut Set", *Electronics* 13 (2024) 603
(`\cite{yeh2024new}`), Figure 7 (test networks 1–22) and Table 9.

All systems are two-terminal networks: the components are the edges, and the
system works iff the source (`Start`) and the target (`End`) are connected by
working edges. The number of components `n` of a system is its number of
**edges**, not of nodes.

## Files

- `yeh2024_raw.txt`: `SistemasYeh_corrected.txt`, the edge lists of the 22
  networks used for this paper since its first version. Kept unchanged, for
  provenance.
- `yeh2024_systems.csv`: one row per system:
  `system, name, label, directed, source, target, n_components, n_nodes`.
- `yeh2024_edges.csv`: one row per component: `system, component, from, to`.
  Component `e` of a system is its `e`-th edge in the raw file; in a state,
  bit `e-1` is 1 iff component `e` works. Do not reorder.

## Reading of the raw file used in this paper

All 22 networks are kept, and they are read exactly as in the raw file:

1. **Networks 6D and 7D are directed**: their edges are arcs, traversable
   only from the first to the second endpoint listed. All other networks are
   undirected.
   - 6D has 4 **irrelevant** components (10–13: the arcs leaving node `4` of
     the raw file, which no path from `Start` reaches), so it is
     semi-coherent but not coherent. Every component of every other system is
     relevant: checked by exhaustive enumeration in the tests for n ≤ 21, and,
     for all undirected systems, by checking that every edge lies in the same
     biconnected component as an added `Start`–`End` edge (i.e., on a simple
     `Start`–`End` path).
   - 6D has 13 minimal cuts and 7D has 110; read as undirected they would
     have 26 and 133 (Table 9 of Yeh et al.: 25 and 133).
2. **Network 21 is network 13 redrawn**: the two are isomorphic with the
   terminals fixed (Section 6.3 of Yeh et al.), so they have the same
   signature and every result for 21 equals the one for 13.
3. **Network 22 is network 18 plus one edge** (`{"2", "3"}` in the raw
   numbering of 22): 23 components instead of 22. Figure 7(22) of Yeh et al.
   indeed has 23 edges, although their Table 9 reports 22 arcs for it.

The companion paper (repository `Simple-delayed`) reads the same file
differently: it drops networks 21 and 22 and treats 6 and 7 as undirected.
Networks 1–5 and 8–20 are identical in both repositories (same edge lists,
same numbering, undirected).

## Check against Table 9 of Yeh et al. (2024)

`n` = nodes, `m` = edges (components), `c` = number of minimal cuts.

| system | ours n, m, c | Table 9 n, m, c | |
|---|---|---|---|
| 1 | 4, 5, 4 | 4, 5, 4 | ok |
| 2 | 6, 8, 9 | 6, 8, 9 | ok |
| 3 | 5, 8, 8 | 5, 8, 8 | ok |
| 4 | 6, 9, 9 | 6, 9, 9 | ok |
| 5 | 9, 12, 28 | 9, 12, 28 | ok |
| 6D | 7, 14, 13 | 7, 14, 25 | directed here (undirected: 26) |
| 7D | 11, 21, 110 | 11, 21, 133 | directed here (undirected: 133) |
| 8 | 9, 13, 24 | 9, 13, 24 | ok |
| 9 | 8, 12, 19 | 8, 12, 19 | ok |
| 10 | 8, 12, 20 | 9, 14, 25 | Table 9 row repeats row 15 |
| 11 | 7, 12, 20 | 7, 12, 20 | ok |
| 12 | 8, 13, 29 | 8, 13, 29 | ok |
| 13 | 16, 30, 644 | 16, 30, 644 (row 21) | ok |
| 14 | 21, 26, 528 | 21, 26, 528 | ok |
| 15 | 9, 14, 25 | 9, 14, 25 | ok |
| 16 | 10, 21, 104 | 10, 21, 104 | ok |
| 17 | 18, 27, 1249 | 18, 27, 1249 | ok |
| 18 | 13, 22, 219 | 13, 22, 214 (row 22) | see Simple-delayed |
| 19 | 20, 30, 7376 | 20, 30, 7376 | ok |
| 20 | 16, 24, 105 | 16, 24, 105 | ok |
| 21 | 16, 30, 644 | = 13 | isomorphic to 13 |
| 22 | 13, 23, — | 13, 22, 214 | 18 plus one edge |

Minimal cuts for systems with m ≥ 22 come from the `Simple-delayed`
repository (C++ enumeration); the Julia tests check those with m ≤ 16.
