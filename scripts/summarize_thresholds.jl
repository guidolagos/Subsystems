# Collect output/thresholds/row_<id>.csv into output/thresholds_summary.csv.
#
# Usage (from the repository root; light): julia --project=. scripts/summarize_thresholds.jl
# Reports missing rows; refuses files of another thresholds_version.

using Subsystems

function main()
    grid = read_grid()
    systems = load_systems()
    out = joinpath(Subsystems.OUTPUT_DIR, "thresholds_summary.csv")
    io = IOBuffer()
    println(io, "row_id,system,name,n,model,beta,mttf_new,n_classes,", join(THRESH_HEADER, ","))
    missing_rows = Int[]
    for r in grid
        path = thresh_path(r.row_id)
        if !isfile(path)
            push!(missing_rows, r.row_id)
            continue
        end
        meta, header, rows = read_threshold_output(path)
        header == THRESH_HEADER || error("$(path): unexpected header")
        (parse(Int, meta["system"]), meta["model"], meta["beta"]) == (r.system, r.model, r.beta) ||
            error("$(path): metadata does not match the grid row")
        sys = get_system(systems, r.system)
        for f in rows
            println(io, r.row_id, ",", r.system, ",", sys.name, ",", meta["n"], ",", r.model, ",", r.beta, ",",
                    meta["mttf_new"], ",", meta["n_classes"], ",", join(f, ","))
        end
    end
    Subsystems.atomic_write(out, take!(io))
    println("wrote $(out): $(length(grid) - length(missing_rows)) of $(length(grid)) rows")
    isempty(missing_rows) || println("missing rows: ", missing_rows)
end

main()
