# Compute structural signatures of two-terminal systems by exhaustive enumeration.
#
# Usage (from the repository root):
#   julia --project=. scripts/compute_signatures.jl [options] [ids...]
#
#   ids        system ids to process, e.g. `7 13 19`; `all` or nothing = all systems
#   --family=F system family (default yeh2024): reads data/systems/F_systems.csv
#   --force    recompute even if the output file already exists
#
# Output: data/signatures/<name>.csv, one file per system, written atomically.
# Existing files are skipped, so the script can be re-run after an interruption.
# Threads: start Julia with -t N (or JULIA_NUM_THREADS=N) to use N threads.

using Subsystems

function main(args)
    force = "--force" in args
    family = "yeh2024"
    ids = Int[]
    for a in args
        if a == "--force"
            continue
        elseif startswith(a, "--family=")
            family = a[length("--family=")+1:end]
        elseif startswith(a, "--")
            error("unknown option $(a)")
        elseif a != "all"
            push!(ids, parse(Int, a))
        end
    end

    systems = load_systems(family)
    selected = isempty(ids) ? systems : [get_system(systems, id) for id in ids]
    println("family = $(family); systems = $([s.id for s in selected]); threads = $(Threads.nthreads())")
    flush(stdout)

    for sys in selected
        path = signature_path(sys)
        if isfile(path) && !force
            println("[skip] $(sys.name): $(path) already exists")
            flush(stdout)
            continue
        end
        n = ncomponents(sys)
        println("[run ] $(sys.name) (label $(sys.label), n = $(n), 2^n = $(2^n) states, ",
                sys.directed ? "directed" : "undirected", ")")
        flush(stdout)
        cs = compile_system(sys)
        t0 = time()
        W = count_working_states(cs)
        seconds = time() - t0
        write_signature_file(path, sys, W; seconds = seconds, threads = Threads.nthreads())
        s = signature(W)
        println("[done] $(sys.name): $(round(seconds; digits = 2)) s; W = $(W)")
        println("       s = ", join((round(Float64(x); digits = 6) for x in s), ", "))
        println("       wrote $(path)")
        flush(stdout)
    end
end

main(ARGS)
