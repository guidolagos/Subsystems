# Benchmark policies (src/policies.jl) for the rows of data/grid.csv.
#
# Usage (from the repository root):
#   julia --project=. scripts/run_policies.jl <task> <ntasks> [options]
#
#   --min-n=a, --max-n=b   only the rows of systems with a <= n <= b (default: all)
#   --rows=a:b             only these row ids
#   --state-max-n=m        optimal state-based policy for n <= m (default 23)
#   --out=main|large       output/policies (default) or output/policies_large
#   --verbose              print the passes of every cost pair
#
# The selected rows are dealt round-robin to the `ntasks` tasks; this task does
# the ones of number `task` (1-based). Rows whose output exists are skipped.
# Exits with code 1 if any row failed.

using Subsystems

function main(args)
    length(args) >= 2 || error("usage: run_policies.jl <task> <ntasks> [options]")
    task, ntasks = parse(Int, args[1]), parse(Int, args[2])
    1 <= task <= ntasks || error("task must be in 1..ntasks")
    minn, maxn, smax, out, verbose, rows = 0, typemax(Int), 23, POLICIES_DIR, false, nothing
    for a in args[3:end]
        if startswith(a, "--min-n=")
            minn = parse(Int, a[9:end])
        elseif startswith(a, "--max-n=")
            maxn = parse(Int, a[9:end])
        elseif startswith(a, "--state-max-n=")
            smax = parse(Int, a[15:end])
        elseif a == "--out=main"
            out = POLICIES_DIR
        elseif a == "--out=large"
            out = POLICIES_LARGE_DIR
        elseif startswith(a, "--rows=")
            lo, hi = parse.(Int, split(a[8:end], ':'))
            rows = lo:hi
        elseif a == "--verbose"
            verbose = true
        else
            error("unknown argument $(a)")
        end
    end
    systems = load_systems()
    nof = Dict(s.id => ncomponents(s) for s in systems)
    sel = [r for r in read_grid() if minn <= nof[r.system] <= maxn && (rows === nothing || r.row_id in rows)]
    mine = [r for (i, r) in enumerate(sel) if (i - 1) % ntasks + 1 == task]
    println("task $(task)/$(ntasks): $(length(mine)) rows; state-based for n <= $(smax); output $(out)")
    flush(stdout)
    done = skipped = failed = 0
    cache = Dict()
    current = 0
    for r in mine
        if r.system != current
            empty!(cache)            # one system in memory at a time
            current = r.system
        end
        t0 = time()
        try
            st = run_policies_row(r; systems = systems, dir = out, state_max_n = smax, cache = cache, verbose = verbose)
            if st == :skipped
                skipped += 1
            else
                done += 1
                println("[row $(r.row_id)] system $(r.system) $(r.model) β=$(r.beta): ", round(time() - t0; digits = 1), " s")
            end
        catch err
            failed += 1
            println("[row $(r.row_id)] FAILED: ", sprint(showerror, err))
        end
        flush(stdout)
    end
    println("done $(done), skipped $(skipped), failed $(failed)")
    failed == 0 || exit(1)
end

main(ARGS)
