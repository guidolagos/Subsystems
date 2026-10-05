# Instructions for coding agents

- Language: Julia 1.10. Activate the project with `--project=.`.
- Put reusable logic in `src/`; keep `scripts/` as thin wrappers.
- Run tests with `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Never commit anything in `output/` or `logs/`.
- Never run heavy computations on the cluster login node; use `sbatch` (or `srun` for setup and tests).
- Everything is exact computation; if randomness is ever needed, use `StableRNGs` seeded from the grid row being processed, never from the Slurm task ID.
- Each unit of work writes one output file atomically (write to a temp file, then rename) and is skipped if that file already exists, so jobs are resumable.
- Systems are two-terminal networks whose components are the edges (`data/systems/<family>_edges.csv`); component `e` is bit `e-1` of a state, and bit = 1 means working. Do not reorder edges: signature files and later experiments rely on this numbering. Systems 6D and 7D are directed.
- Signatures of the complete systems live in `data/signatures/<name>.csv` and are committed. They are computed by `scripts/compute_signatures.jl` (cluster: `templates/signatures.sbatch`) and validated by `scripts/check_signatures.jl` against `test/reference/`.
- Output files carry a version number; bump it whenever their content changes, and make readers refuse other versions.
- `test/reference/*.py` are independent Python implementations used only to generate the reference files compared by `Pkg.test()`; they are not run on the cluster.
