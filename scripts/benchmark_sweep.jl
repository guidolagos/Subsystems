# Computational effort of Algorithm alg (per-state signatures and MTTF) against the
# enumeration of the 2^n states needed for the signature of the complete system.
#
# Usage (from the repository root):  julia --project=. scripts/benchmark_sweep.jl <system id>...
#
# For each system, on one thread, after a warm-up on system 1:
#   seconds_enum      count_working_states: Φ at the 2^n states, counts W_k (signature of the system);
#   seconds_sweep     Algorithm alg: counting vectors W^x of all 2^n states, grouped into classes;
#   seconds_sweep_mttf  the same sweep, also computing at every state the MTTF under the 31
#                     failure models of the grid (Exp and 6 models x 5 betas);
#   seconds_tables    computing the 31 model tables (δ^w_k, τ_w).
# The per-state MTTF sums are checked against the class-based computation.
# Output: output/benchmark/sweep_<name>.csv. Peak memory: see sacct MaxRSS.

using Subsystems

function grid_tables(n)
    ms = [lfmo_model("Exp", "0")]
    for name in MODELS[2:end], b in GRID_BETAS
        push!(ms, lfmo_model(name, b))
    end
    return [model_tables(m, n) for m in ms]
end

function bench(sys)
    cs = compile_system(sys)
    n = cs.n
    t = time(); W = count_working_states(cs; nchunks = 1); t_enum = time() - t
    t = time(); tabs = grid_tables(n); t_tab = time() - t
    res = sweep(cs)
    check_sweep(res, W)
    res2 = sweep(cs; deltas = [mt.delta for mt in tabs])
    # per-state sums vs classes
    for (i, mt) in enumerate(tabs)
        s = 0.0
        for w in 1:n, (c, k) in enumerate(res.levels[w + 1].mult)
            s += k * class_mttf(mt, view(res.levels[w + 1].W, :, c), w)
        end
        isapprox(s, res2.mttf_sums[i]; rtol = 1e-9) || error("$(sys.name): per-state MTTF sum differs from classes (model $(i))")
    end
    return (n = n, t_enum = t_enum, t_sweep = res.seconds, t_sweep_mttf = res2.seconds, t_tab = t_tab,
            working = sum(res.working), classes = nclasses(res), mps = res.mps, bytes = res.level_bytes,
            evals = res.phi_evals)
end

function main(args)
    systems = load_systems()
    ids = isempty(args) ? [s.id for s in systems] : parse.(Int, args)
    bench(get_system(systems, 1))                          # warm-up (compilation)
    mkpath(Subsystems.BENCH_DIR)
    for id in ids
        sys = get_system(systems, id)
        println("[bench] $(sys.name) (n = $(ncomponents(sys)))"); flush(stdout)
        b = bench(sys)
        io = IOBuffer()
        println(io, "name,n,seconds_enum,seconds_sweep,seconds_sweep_mttf,seconds_tables,phi_evals,",
                "working_states,classes,mps,level_bytes,julia_version,host")
        println(io, sys.name, ",", b.n, ",", round(b.t_enum; digits = 3), ",", round(b.t_sweep; digits = 3), ",",
                round(b.t_sweep_mttf; digits = 3), ",", round(b.t_tab; digits = 3), ",", b.evals, ",", b.working, ",",
                b.classes, ",", b.mps, ",", b.bytes, ",", VERSION, ",", Subsystems._hostname())
        path = joinpath(Subsystems.BENCH_DIR, "sweep_$(sys.name).csv")
        Subsystems.atomic_write(path, take!(io))
        println("        enum $(round(b.t_enum; digits = 2)) s, sweep $(round(b.t_sweep; digits = 2)) s, ",
                "sweep + 31 MTTF $(round(b.t_sweep_mttf; digits = 2)) s, $(b.classes) classes; wrote $(path)")
        flush(stdout)
    end
end

main(ARGS)
