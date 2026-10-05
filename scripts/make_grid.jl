# Write data/grid.csv: one row per (system, failure model, beta); 22 x 31 = 682 rows.
#
# Usage (from the repository root): julia --project=. scripts/make_grid.jl

using Subsystems

rows = make_grid()
write_grid(GRID_PATH, rows)
println("wrote $(GRID_PATH): $(length(rows)) rows ($(ROWS_PER_SYSTEM) per system)")
