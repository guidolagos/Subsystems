# MTTF-threshold policies for the rows of one system (src/grid.jl).
#
# Usage (from the repository root):
#   julia --project=. scripts/run_system.jl <system id> [--rows=a:b]
#
# The states are swept once (cached in output/classes/<name>.bin); then every pending
# row of the system writes output/thresholds/row_<id>*.csv. Rows whose main file
# exists are skipped, so the script can be re-run after an interruption.
# Exits with code 1 if any row failed.

using Subsystems

function main(args)
    isempty(args) && error("usage: run_system.jl <system id> [--rows=a:b]")
    sysid = parse(Int, args[1])
    rows = nothing
    for a in args[2:end]
        if startswith(a, "--rows=")
            lo, hi = parse.(Int, split(a[length("--rows=")+1:end], ':'))
            rows = lo:hi
        else
            error("unknown argument $(a)")
        end
    end
    t0 = time()
    done, skipped, failed = run_system(sysid; rows = rows)
    println("system $(sysid): done $(done), skipped $(skipped), failed $(failed); ",
            round(time() - t0; digits = 1), " s")
    failed == 0 || exit(1)
end

main(ARGS)
