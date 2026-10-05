using Test
using Subsystems

# Build a small system from an edge list.
mksys(edges; directed = false, name = "test") =
    TwoTerminalSystem(0, name, name, directed, "Start", "End", [(String(u), String(v)) for (u, v) in edges])

sig(edges; directed = false) = signature(count_working_states(mksys(edges; directed = directed)))

# Independent structure function: plain depth-first search on adjacency lists.
function naive_works(sys::TwoTerminalSystem, x::UInt64)
    adj = Dict{String,Vector{String}}()
    for (e, (u, v)) in enumerate(sys.edges)
        (x >> (e - 1)) & 1 == 1 || continue
        push!(get!(adj, u, String[]), v)
        sys.directed || push!(get!(adj, v, String[]), u)
    end
    seen = Set([sys.source])
    stack = [sys.source]
    while !isempty(stack)
        u = pop!(stack)
        u == sys.target && return true
        for v in get(adj, u, String[])
            v in seen || (push!(seen, v); push!(stack, v))
        end
    end
    return false
end

# Deterministic pseudo-random 64-bit states (xorshift64), no RNG dependency.
function xorshift_states(seed::UInt64, count::Integer, n::Integer)
    x = seed
    out = UInt64[]
    mask = (UInt64(1) << n) - one(UInt64)
    for _ in 1:count
        x ⊻= x << 13; x ⊻= x >> 7; x ⊻= x << 17
        push!(out, x & mask)
    end
    return out
end

# Number of minimal cut sets: failed sets F such that the system fails, and
# repairing any single component of F makes it work. Exhaustive; small n only.
function count_minimal_cuts(sys::TwoTerminalSystem)
    n = ncomponents(sys)
    cs = compile_system(sys)
    N = UInt64(1) << n
    phi = [works(cs, x) for x in UInt64(0):(N - one(UInt64))]
    full = N - one(UInt64)
    count = 0
    for x in UInt64(0):full
        phi[x + 1] && continue
        failed = full & ~x
        minimal = true
        for e in 0:(n - 1)
            if (failed >> e) & 1 == 1 && !phi[(x | (UInt64(1) << e)) + 1]
                minimal = false
                break
            end
        end
        count += minimal
    end
    return count
end

function read_reference(path)
    ref = Dict{Int,Vector{Int}}()
    for line in eachline(path)
        (startswith(line, "#") || startswith(line, "system") || isempty(strip(line))) && continue
        f = split(line, ',')
        push!(get!(ref, parse(Int, f[1]), Int[]), parse(Int, f[3]))
    end
    return ref
end

# Components whose state matters for at least one state of the others (exhaustive).
function relevant_components(sys::TwoTerminalSystem)
    n = ncomponents(sys)
    cs = compile_system(sys)
    N = UInt64(1) << n
    phi = [works(cs, x) for x in UInt64(0):(N - one(UInt64))]
    rel = falses(n)
    for x in UInt64(0):(N - one(UInt64)), e in 1:n
        b = UInt64(1) << (e - 1)
        x & b == zero(UInt64) || continue
        phi[x + 1] != phi[(x | b) + 1] && (rel[e] = true)
    end
    return findall(rel)
end

# Naive counting vector of a state: enumerate its sub-states directly.
function naive_counts(cs::CompiledSystem, x::UInt64)
    w = count_ones(x)
    W = zeros(Int, w + 1)
    y = x
    while true
        works(cs, y) && (W[w - count_ones(y) + 1] += 1)
        y == 0 && break
        y = (y - one(UInt64)) & x
    end
    return W
end

function read_thresholds_reference(path)
    header, rows = Subsystems.read_simple_csv(path)
    col(name) = findfirst(==(name), header)
    cases = Dict{Tuple{Int,String,String},Vector{Vector{String}}}()
    for r in rows
        push!(get!(cases, (parse(Int, r[col("system")]), r[col("model")], r[col("beta")]), Vector{String}[]), r)
    end
    return header, cases
end

@testset "Subsystems" begin

    @testset "known signatures" begin
        # series system of 4 components: fails at the first failure
        @test sig([("Start", "a"), ("a", "b"), ("b", "c"), ("c", "End")]) == [1, 0, 0, 0]
        # parallel system of 3 components: fails at the last failure
        @test sig([("Start", "End"), ("Start", "End"), ("Start", "End")]) == [0, 0, 1]
        # component 1 in series with the parallel pair (2, 3)
        @test sig([("Start", "a"), ("a", "End"), ("a", "End")]) == [1 // 3, 2 // 3, 0]
        # bridge system (Samaniego 2007): s = (0, 1/5, 3/5, 1/5, 0)
        bridge = [("Start", "1"), ("Start", "2"), ("1", "2"), ("1", "End"), ("2", "End")]
        @test sig(bridge) == [0, 1 // 5, 3 // 5, 1 // 5, 0]
        # directedness is honoured: an arc End -> Start is irrelevant when directed
        @test sig([("Start", "End"), ("End", "Start")]; directed = true) == [1 // 2, 1 // 2]
        @test sig([("Start", "End"), ("End", "Start")]; directed = false) == [0, 1]
        # dangling component (semi-coherent, irrelevant component)
        @test sig([("Start", "a"), ("a", "End"), ("b", "a")]) == [2 // 3, 1 // 3, 0]
    end


    @testset "state convention and counts" begin
        bridge = mksys([("Start", "1"), ("Start", "2"), ("1", "2"), ("1", "End"), ("2", "End")])
        @test works(bridge, state_from_working(5, [1, 4]))          # path Start-1-End
        @test works(bridge, state_from_working(5, [2, 3, 4]))       # path Start-2-1-End
        @test !works(bridge, state_from_working(5, [1, 2, 3]))      # nothing reaches End
        @test count_working_states(bridge) == [1, 5, 8, 2, 0, 0]
        # result independent of the chunking of the state space
        @test count_working_states(bridge; nchunks = 1) == count_working_states(bridge; nchunks = 7)
        @test_throws ErrorException signature([0, 0, 0])            # Φ(1) = 0 is rejected
    end

    systems = load_systems("yeh2024")

    @testset "Yeh et al. (2024) data files" begin
        @test [s.id for s in systems] == collect(1:22)
        @test [ncomponents(s) for s in systems] ==
              [5, 8, 8, 9, 12, 14, 21, 13, 12, 12, 12, 13, 30, 26, 14, 21, 27, 22, 30, 24, 30, 23]
        @test [nnodes(s) for s in systems] ==
              [4, 6, 5, 6, 9, 7, 11, 9, 8, 8, 7, 8, 16, 21, 9, 10, 18, 13, 20, 16, 16, 13]
        # systems 6D and 7D are directed (arcs as listed in the raw file); the others undirected
        @test [s.id for s in systems if s.directed] == [6, 7]
        @test [s.label for s in systems if s.directed] == ["6D", "7D"]
        # the CSV edge lists coincide with the raw file, for all 22 systems
        raw = parse_yeh_raw(joinpath(Subsystems.SYSTEMS_DIR, "yeh2024_raw.txt"))
        @test [r.id for r in raw] == collect(1:22)
        for s in systems
            r = raw[s.id]
            @test (r.id, r.label, r.directed, r.edges) == (s.id, s.label, s.directed, s.edges)
        end
    end

    @testset "minimal cuts, systems with n <= 16 (see data/systems/README.md)" begin
        # Table 9 of Yeh et al. (2024), except system 6 (directed here: 13; Table 9: 25,
        # undirected: 26) and system 10 (Table 9 repeats the row of system 15)
        expected = Dict(1 => 4, 2 => 9, 3 => 8, 4 => 9, 5 => 28, 6 => 13, 8 => 24,
                        9 => 19, 10 => 20, 11 => 20, 12 => 29, 15 => 25)
        for sys in systems
            ncomponents(sys) <= 16 || continue
            @test count_minimal_cuts(sys) == expected[sys.id]
        end
    end

    @testset "relevant components, systems with n <= 21" begin
        # only system 6D is semi-coherent: its components 10-13 (the arcs leaving node 4
        # of the raw file, which no path from Start reaches) are irrelevant
        for sys in systems
            ncomponents(sys) <= 21 || continue
            expected = sys.id == 6 ? setdiff(1:14, 10:13) : collect(1:ncomponents(sys))
            @test relevant_components(sys) == expected
        end
    end

    @testset "structure function vs. naive DFS (random states, all Yeh systems)" begin
        for sys in systems
            cs = compile_system(sys)
            states = xorshift_states(UInt64(0x9E3779B97F4A7C15) + UInt64(sys.id), 3000, ncomponents(sys))
            @test all(works(cs, x) == naive_works(sys, x) for x in states)
        end
    end


    @testset "counts vs. independent reference (Yeh systems with n <= 23)" begin
        ref = read_reference(joinpath(@__DIR__, "reference", "yeh2024_W_reference.csv"))
        @test sort(collect(keys(ref))) == collect(1:22)
        @test ref[21] == ref[13]                       # isomorphic networks
        for sys in systems
            @test length(ref[sys.id]) == ncomponents(sys) + 1
            ncomponents(sys) <= 23 || continue
            W = count_working_states(sys)
            @test W == ref[sys.id]
            @test sum(signature(W)) == 1
        end
    end

    @testset "signature file round trip" begin
        sys = get_system(systems, 1)
        W = count_working_states(sys)
        dir = mktempdir()
        path = joinpath(dir, "bridge.csv")
        write_signature_file(path, sys, W; seconds = 0.0, threads = 1)
        r = read_signature_file(path)
        @test r.W == W
        @test r.s == [0, 1 // 5, 3 // 5, 1 // 5, 0]
        @test r.meta["name"] == "yeh2024_01"
        @test readdir(dir) == ["bridge.csv"]          # no leftover temporary file
    end


    @testset "LFMO models: calibration by the pairwise beta factor" begin
        for name in MODELS[2:end], b in ("0.01", "0.05", "0.1", "0.2", "0.5")
            m = lfmo_model(name, b)
            @test abs(psi(m, 1) - 1) < 1e-60
            @test abs(2 - psi(m, 2) - parse(BigFloat, b)) < 1e-60
            @test psi(m, 0) == 0
        end
        # closed forms of the parameter
        @test Float64(lfmo_model("LFMO-C", "0.2").param) ≈ 0.2 / 1.8
        @test Float64(lfmo_model("BFR", "0.2").param) == 0.2
        @test Float64(lfmo_model("LSh", "0.2").param) == 0.2
        @test Float64(lfmo_model("LFMO-S", "0.5").param) ≈ log2(1.5)
        # beta = 1/2 gives eta = (sqrt(5)-1)/2 for Gamma and 1/sqrt(12) for inverse Gaussian
        @test Float64(lfmo_model("LFMO-G", "0.5").param) ≈ (sqrt(5) - 1) / 2
        @test Float64(lfmo_model("LFMO-IG", "0.5").param) ≈ 1 / sqrt(12)
        # the earlier parameters of the paper: Ψ(x) = 10x/(x+9) is LFMO-C with p = 0.1
        @test Float64(lfmo_model("LFMO-C", string(2 - 20 / 11)).param) ≈ 0.1
        # inverse Gaussian cannot exceed beta = 2 - sqrt(2)
        @test_throws ErrorException lfmo_model("LFMO-IG", "0.6")
        @test_throws ErrorException lfmo_model("Exp", "0.1")
        @test_throws ErrorException lfmo_model("LFMO-C", "0")
    end


    @testset "transition rates of the survivor-count chain" begin
        n = 30
        for name in MODELS, b in (name == "Exp" ? ("0",) : ("0.01", "0.5"))
            m = lfmo_model(name, b)
            R, ex = transition_rates(m, n)
            @test all(>=(0), R)
            for c in 0:(n - 1)
                @test isapprox(sum(R[c + 1, :]), ex[c + 1]; rtol = 1e-13)
                @test ex[c + 1] ≈ Float64(psi(m, n - c))
                # mean batch = w Ψ(1) / Ψ(w)
                w = n - c
                @test isapprox(sum((f) * R[c + 1, c + f + 1] for f in 1:w) / ex[c + 1], mean_batch(m, w); rtol = 1e-12)
            end
        end
        # iid exponential: single failures only
        R, ex = transition_rates(lfmo_model("Exp", "0"), 6)
        @test all(R[c + 1, c2 + 1] == (c2 == c + 1 ? 6 - c : 0) for c in 0:5, c2 in 0:6)
        # lethal shock: drift (1-θ) per component plus killing θ of all survivors
        R, ex = transition_rates(lfmo_model("LSh", "0.2"), 5)
        @test R[1, 2] ≈ 5 * 0.8 && R[1, 6] ≈ 0.2 && all(abs(R[1, j]) < 1e-15 for j in 3:5)
        # BFR: binomial batches, shock rate 1/p
        R, ex = transition_rates(lfmo_model("BFR", "0.2"), 4)
        @test all(R[1, f + 1] ≈ binomial(4, f) * 0.2^f * 0.8^(4 - f) / 0.2 for f in 1:4)
    end


    @testset "sweep: counting vectors of all states (Lemma lemma:Wk, Algorithm alg)" begin
        for id in (1, 3, 6, 9)
            sys = get_system(systems, id)
            cs = compile_system(sys)
            n = cs.n
            rec = Dict{UInt64,Vector{Int}}()
            res = sweep(cs; record = rec)
            @test length(rec) == 2^n
            @test all(rec[x] == naive_counts(cs, x) for x in UInt64(0):UInt64(2^n - 1))
            @test res.phi_evals == 2^n
            @test check_sweep(res, count_working_states(cs))
            # classes = distinct counting vectors of working states, with multiplicities
            mult = Dict{Vector{Int},Int}()
            for (x, W) in rec
                works(cs, x) && (mult[W] = get(mult, W, 0) + 1)
            end
            @test nclasses(res) == length(mult)
            @test all(mult[Int.(res.levels[w + 1].W[:, c])] == res.levels[w + 1].mult[c]
                      for w in 0:n for c in eachindex(res.levels[w + 1].mult))
            # minimal path sets: working states whose sons all fail
            mps = count(x -> works(cs, x) && all(!works(cs, x ⊻ (UInt64(1) << e)) for e in 0:(n - 1) if (x >> e) & 1 == 1),
                        UInt64(0):UInt64(2^n - 1))
            @test res.mps == mps
        end
        # known cases
        bridge = mksys([("Start", "1"), ("Start", "2"), ("1", "2"), ("1", "End"), ("2", "End")])
        @test sweep(bridge).mps == 4
        @test sweep(mksys([("Start", "a"), ("a", "b"), ("b", "End")])).mps == 1
        @test sweep(mksys([("Start", "End"), ("Start", "End"), ("Start", "End")])).mps == 3
        r = sweep(mksys([("Start", "a"), ("a", "b"), ("b", "End")]))
        @test nclasses(r) == 1 && Int.(r.levels[4].W[:, 1]) == [1, 0, 0, 0]     # series: only (1,1,1) works
    end

    @testset "sweep: systems with n <= 23 agree with the committed signatures; class files" begin
        for sys in systems
            ncomponents(sys) <= 23 || continue
            res = sweep(sys)
            @test check_sweep(res, read_signature_file(signature_path(sys)).W)
            @test sum(res.working) == sum(read_signature_file(signature_path(sys)).W)
        end
        res = sweep(get_system(systems, 6))
        dir = mktempdir()
        path = joinpath(dir, "c.bin")
        write_classes(path, res)
        r2 = read_classes(path)
        @test r2.n == res.n && r2.mps == res.mps && r2.working == res.working && r2.phi_evals == res.phi_evals
        @test all(r2.levels[w].W == res.levels[w].W && r2.levels[w].mult == res.levels[w].mult for w in eachindex(res.levels))
        @test readdir(dir) == ["c.bin"]
    end

    @testset "model tables: δ, τ and order statistics" begin
        n = 12
        mt = model_tables(lfmo_model("Exp", "0"), n)
        @test all(isapprox(mt.delta[w, k + 1], 1 / (w - k); rtol = 1e-14) for w in 1:n for k in 0:(w - 1))
        @test all(isapprox(mt.tau[w], 1 / (w * binomial(n, w)); rtol = 1e-13) for w in 1:n)
        for name in MODELS, b in (name == "Exp" ? ("0",) : ("0.01", "0.5"))
            m = lfmo_model(name, b)
            mt = model_tables(m, n)
            R, ex = transition_rates(m, n)
            # E T_{k:w} by the recursion eq:ETkn of the paper (nonnegative terms)
            ET = zeros(n, n)
            for w in 1:n, k in 1:w
                s = 1 / ex[n - w + 1]
                for j in 1:(k - 1)
                    s += R[n - w + 1, n - w + j + 1] / ex[n - w + 1] * ET[w - j, k - j]
                end
                ET[w, k] = s
            end
            @test all(isapprox(expected_order_statistics(mt, w), ET[w, 1:w]; rtol = 1e-11) for w in 1:n)
            @test mt.tau[n] ≈ 1 / Float64(psi(m, n))
            @test all(mt.q[w, j] >= 0 for w in 1:n for j in 1:w)
            @test all(isapprox(sum(binomial(w, j) * mt.q[w, j] for j in 1:w), mt.psi[w]; rtol = 1e-12) for w in 1:n)
        end
    end

    @testset "threshold curves vs independent reference (test/reference/thresholds_reference.csv)" begin
        header, cases = read_thresholds_reference(joinpath(@__DIR__, "reference", "thresholds_reference.csv"))
        col(name) = findfirst(==(name), header)
        cache = Dict{Int,SweepResult}()
        @test length(cases) == 24
        for ((id, model, beta), rows) in cases
            sys = get_system(systems, id)
            res = get!(() -> sweep(sys), cache, id)
            cv = threshold_curve(res, model_tables(lfmo_model(model, beta), res.n))
            @test length(cv.gamma) == parse(Int, rows[1][col("n_breakpoints")])
            @test res.mps == parse(Int, rows[1][col("mps")])
            @test isapprox(cv.m1, parse(Float64, rows[1][col("m1")]); rtol = 1e-11)
            for r in rows
                i = parse(Int, r[col("index")])
                for (got, name) in ((cv.gamma[i], "gamma"), (cv.T[i], "T"), (cv.N[i], "N"), (cv.P[i], "P"))
                    @test isapprox(got, parse(Float64, r[col(name)]); rtol = 1e-10, atol = 1e-14)
                end
            end
        end
    end

    @testset "threshold curves: identities, costs, local minima" begin
        sys = get_system(systems, 1)                  # bridge
        res = sweep(sys)
        cv = threshold_curve(res, model_tables(lfmo_model("Exp", "0"), 5))
        @test cv.m1 ≈ 49 / 60                          # MTTF of the bridge with iid Exp(1)
        @test cv.T[1] ≈ 1 / 5 && cv.N[1] ≈ 1 && cv.P[1] == 0     # γ = 1: repair at the first failure
        @test cv.T[end] ≈ cv.m1 && cv.P[end] ≈ 1 && cv.states[end] == 16
        L = curve_ltmc(cv; ratio = 2.0, setup = 0.1)
        @test L[1] ≈ (0.1 + 2.0 / 5 * 1) / (1 / 5)
        @test curve_at(cv, 0) == length(cv.gamma) && curve_at(cv, 1) == 1
        @test curve_at(cv, cv.gamma[3]) == 3 && curve_at(cv, cv.gamma[3] + 1e-9) == 2
        @test local_minima([3, 2, 2, 1, 2]) == 1
        @test local_minima([1, 2, 1]) == 2
        @test local_minima([1, 1, 1]) == 1
        # all models on a semi-coherent system: monotone curves of T and P
        res6 = sweep(get_system(systems, 6))
        for name in MODELS, b in (name == "Exp" ? ("0",) : ("0.05", "0.5"))
            cv6 = threshold_curve(res6, model_tables(lfmo_model(name, b), 14))
            @test issorted(cv6.T) && issorted(cv6.states) && all(diff(cv6.gamma) .< 0)
            @test all(0 .<= cv6.P .<= 1 + 1e-12)
        end
    end

    @testset "experiment grid and row outputs" begin
        grid = read_grid()
        @test length(grid) == 22 * ROWS_PER_SYSTEM == 682
        @test grid == make_grid()
        @test all(grid[r].system == (r - 1) ÷ ROWS_PER_SYSTEM + 1 for r in eachindex(grid))
        dir, cdir, sdir = mktempdir(), mktempdir(), mktempdir()
        done, skipped, failed = run_system(1; grid = grid, systems = systems, dir = dir, classes_dir = cdir,
                                           sweep_dir = sdir, rows = 1:3, verbose = false)
        @test (done, skipped, failed) == (3, 0, 0)
        @test isfile(joinpath(cdir, "yeh2024_01.bin")) && isfile(joinpath(sdir, "yeh2024_01.csv"))
        meta, header, rows = read_threshold_output(thresh_path(2; dir = dir))
        @test header == THRESH_HEADER && length(rows) == length(COST_RATIOS) * length(COST_SETUPS)
        @test meta["model"] == "LFMO-C" && meta["beta"] == "0.01" && meta["n"] == "5"
        meta, header, rows = read_threshold_output(curve_path(1; dir = dir))
        @test length(rows) == length(CURVE_GAMMAS)
        @test isfile(breaks_path(3; dir = dir))
        # the best LTMC is the minimum over the breakpoints, below γ = 1 and corrective
        _, h, rr = read_threshold_output(thresh_path(1; dir = dir))
        for f in rr
            best, g1, corr = parse.(Float64, f[findfirst(==("ltmc_best"), h):findfirst(==("ltmc_corrective"), h)])
            @test best <= g1 && best <= corr
        end
        @test run_system(1; grid = grid, systems = systems, dir = dir, classes_dir = cdir, sweep_dir = sdir,
                         rows = 1:3, verbose = false) == (0, 3, 0)
    end

end
