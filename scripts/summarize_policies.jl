# Collect output/policies/row_<id>.csv (and output/policies_large/, whose state-based
# columns take precedence) into output/policies_summary.csv, one line per row and
# cost pair, with the relative gaps
#   gap_gamma_state = ltmc_best_gamma / ltmc_state - 1   (best MTTF-threshold vs optimal state-based)
#   gap_count_state = ltmc_count / ltmc_state - 1
#   gap_r_gamma     = ltmc_best_r / ltmc_best_gamma - 1   (best r-out-of-n vs best MTTF-threshold)
#
# Usage (from the repository root; light): julia --project=. scripts/summarize_policies.jl

using Subsystems

function main()
    grid = read_grid()
    out = joinpath(Subsystems.OUTPUT_DIR, "policies_summary.csv")
    io = IOBuffer()
    println(io, "row_id,system,model,beta,n,", join(POLICIES_HEADER, ","), ",gap_gamma_state,gap_count_state,gap_r_gamma")
    missing_rows = Int[]
    for r in grid
        path = policies_path(r.row_id)
        if !isfile(path)
            push!(missing_rows, r.row_id)
            continue
        end
        meta, header, rows = read_policies_output(path)
        largep = policies_path(r.row_id; dir = POLICIES_LARGE_DIR)
        if isfile(largep)
            _, _, brows = read_policies_output(largep)
            rows = brows
        end
        col(name) = findfirst(==(name), header)
        for f in rows
            Lg = parse(Float64, f[col("ltmc_best_gamma")]); Ls = parse(Float64, f[col("ltmc_state")])
            Lc = parse(Float64, f[col("ltmc_count")]); Lr = parse(Float64, f[col("ltmc_best_r")])
            println(io, r.row_id, ",", r.system, ",", r.model, ",", r.beta, ",", meta["n"], ",", join(f, ","), ",",
                    repr(Lg / Ls - 1), ",", repr(Lc / Ls - 1), ",", repr(Lr / Lg - 1))
        end
    end
    Subsystems.atomic_write(out, take!(io))
    println("wrote $(out): $(length(grid) - length(missing_rows)) of $(length(grid)) rows")
    isempty(missing_rows) || println("missing rows: ", missing_rows)
end

main()
