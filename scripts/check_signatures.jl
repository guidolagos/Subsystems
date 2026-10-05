# Check the computed signature files and summarize them.
#
# Usage (from the repository root; light, can run anywhere):
#   julia --project=. scripts/check_signatures.jl [--family=yeh2024]
#
# For every system of the family it reports whether data/signatures/<name>.csv
# exists, whether its W counts coincide with the independent reference in
# test/reference/<family>_W_reference.csv (computed independently of the
# Julia code; see test/reference/W_reference.py), and the signature. Exits with code 1 if any file is
# missing, inconsistent, or differs from the reference.

using Subsystems

function read_reference(path)
    ref = Dict{Int,Vector{Int}}()
    isfile(path) || return ref
    for line in eachline(path)
        line = strip(line)
        (isempty(line) || startswith(line, "#") || startswith(line, "system")) && continue
        f = split(line, ',')
        id, k, w = parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3])
        v = get!(ref, id, Int[])
        length(v) == k || error("reference file: k values out of order for system $(id)")
        push!(v, w)
    end
    return ref
end

function main(args)
    family = "yeh2024"
    for a in args
        startswith(a, "--family=") ? (family = a[length("--family=")+1:end]) : error("unknown argument $(a)")
    end
    systems = load_systems(family)
    ref_path = joinpath(Subsystems.REPO_DIR, "test", "reference", "$(family)_W_reference.csv")
    ref = read_reference(ref_path)
    isempty(ref) && println("WARNING: no reference file at $(ref_path)")

    nbad = 0
    println(rpad("system", 12), rpad("label", 6), rpad("n", 4), rpad("dir", 5), rpad("status", 25),
            rpad("seconds", 10), "E[index] = sum_k k s_k")
    for sys in systems
        path = signature_path(sys)
        status, seconds, eidx = "", "", ""
        if !isfile(path)
            status = "MISSING"
            nbad += 1
        else
            try
                r = read_signature_file(path)       # checks internal consistency
                length(r.W) == ncomponents(sys) + 1 || error("wrong length")
                if haskey(ref, sys.id)
                    status = r.W == ref[sys.id] ? "ok (= reference)" : "DIFFERS FROM REFERENCE"
                    r.W == ref[sys.id] || (nbad += 1)
                else
                    status = "ok (no reference)"
                end
                seconds = get(r.meta, "seconds", "")
                eidx = string(round(Float64(sum(k * r.s[k] for k in eachindex(r.s))); digits = 4))
            catch err
                status = "INVALID: " * sprint(showerror, err)
                nbad += 1
            end
        end
        println(rpad(sys.name, 12), rpad(sys.label, 6), rpad(string(ncomponents(sys)), 4),
                rpad(sys.directed ? "yes" : "no", 5), rpad(status, 25), rpad(seconds, 10), eidx)
    end
    println()
    if nbad == 0
        println("All $(length(systems)) signature files present, consistent and equal to the reference.")
    else
        println("$(nbad) problem(s) found.")
        exit(1)
    end
end

main(ARGS)
