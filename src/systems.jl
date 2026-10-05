# Two-terminal network systems: data loading.
#
# A system is a (directed or undirected) multigraph whose edges are the
# components. The system works iff there is a path of working edges from
# `source` to `target`. Components are numbered 1..n in the order in which
# they appear in the edge list; this numbering is the one used everywhere
# else in the package (bit e-1 of a state is component e).

"""
    TwoTerminalSystem

A two-terminal network system. `edges[e] = (from, to)` is component `e`.
If `directed` is true, component `e` can only be traversed from `from` to `to`.
"""
struct TwoTerminalSystem
    id::Int
    name::String
    label::String
    directed::Bool
    source::String
    target::String
    edges::Vector{Tuple{String,String}}
end

ncomponents(sys::TwoTerminalSystem) = length(sys.edges)

function nnodes(sys::TwoTerminalSystem)
    nodes = Set{String}()
    for (u, v) in sys.edges
        push!(nodes, u)
        push!(nodes, v)
    end
    return length(nodes)
end

Base.show(io::IO, sys::TwoTerminalSystem) = print(io,
    "TwoTerminalSystem(", sys.name, ", label=", sys.label,
    ", n=", ncomponents(sys), ", ", sys.directed ? "directed" : "undirected", ")")

# ---------------------------------------------------------------------------
# Minimal CSV reading (no quoting, no commas inside fields). Our data files
# are written by us, so this is enough and avoids a CSV.jl dependency.
# ---------------------------------------------------------------------------

"""
    read_simple_csv(path) -> (header::Vector{String}, rows::Vector{Vector{String}})

Read a comma-separated file without quoted fields. Blank lines and lines
starting with `#` are ignored.
"""
function read_simple_csv(path::AbstractString)
    header = String[]
    rows = Vector{Vector{String}}()
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "#")) && continue
        fields = [String(strip(f)) for f in split(line, ',')]
        if isempty(header)
            header = fields
        else
            length(fields) == length(header) || error(
                "$(path): expected $(length(header)) fields, got $(length(fields)) in line: $(raw)")
            push!(rows, fields)
        end
    end
    return header, rows
end

function column(header::Vector{String}, name::AbstractString)
    i = findfirst(==(name), header)
    i === nothing && error("column '$(name)' not found in header $(header)")
    return i
end

"""
    load_systems(family = "yeh2024"; dir = SYSTEMS_DIR) -> Vector{TwoTerminalSystem}

Load the systems of a family from `<dir>/<family>_systems.csv` and
`<dir>/<family>_edges.csv`, sorted by system id. Checks that the number of
components and nodes of each system match the systems table.
"""
function load_systems(family::AbstractString = "yeh2024"; dir::AbstractString = SYSTEMS_DIR)
    hs, rs = read_simple_csv(joinpath(dir, "$(family)_systems.csv"))
    he, re = read_simple_csv(joinpath(dir, "$(family)_edges.csv"))

    edges_by_system = Dict{Int,Vector{Tuple{Int,String,String}}}()
    ci_sys, ci_comp = column(he, "system"), column(he, "component")
    ci_from, ci_to = column(he, "from"), column(he, "to")
    for r in re
        id = parse(Int, r[ci_sys])
        push!(get!(edges_by_system, id, Tuple{Int,String,String}[]),
              (parse(Int, r[ci_comp]), r[ci_from], r[ci_to]))
    end

    systems = TwoTerminalSystem[]
    for r in rs
        id = parse(Int, r[column(hs, "system")])
        list = sort(get(edges_by_system, id, Tuple{Int,String,String}[]); by = first)
        [c for (c, _, _) in list] == collect(1:length(list)) ||
            error("system $(id): components must be numbered 1..n without gaps")
        sys = TwoTerminalSystem(
            id,
            r[column(hs, "name")],
            r[column(hs, "label")],
            parse(Bool, r[column(hs, "directed")]),
            r[column(hs, "source")],
            r[column(hs, "target")],
            [(u, v) for (_, u, v) in list],
        )
        ncomponents(sys) == parse(Int, r[column(hs, "n_components")]) ||
            error("system $(id): n_components in table does not match edge list")
        nnodes(sys) == parse(Int, r[column(hs, "n_nodes")]) ||
            error("system $(id): n_nodes in table does not match edge list")
        push!(systems, sys)
    end
    sort!(systems; by = s -> s.id)
    length(unique(s.id for s in systems)) == length(systems) || error("duplicate system ids")
    return systems
end

"""
    get_system(systems, id) -> TwoTerminalSystem
"""
function get_system(systems::Vector{TwoTerminalSystem}, id::Integer)
    i = findfirst(s -> s.id == id, systems)
    i === nothing && error("system $(id) not found; available ids: $([s.id for s in systems])")
    return systems[i]
end

"""
    parse_yeh_raw(path) -> Vector{NamedTuple}

Parse the original text file with the systems of Yeh et al. (2024), in the
format

    6D:
    {"Start", "1"},
    ...

A trailing `D` in the label marks a directed system. Returns a vector of
`(id, label, directed, edges)`. Used only to check that the CSV files agree
with the original source.
"""
function parse_yeh_raw(path::AbstractString)
    out = NamedTuple{(:id, :label, :directed, :edges),
                     Tuple{Int,String,Bool,Vector{Tuple{String,String}}}}[]
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        m = match(r"^(\d+)(D?):$", line)
        if m !== nothing
            push!(out, (id = parse(Int, m.captures[1]), label = String(chop(line)),
                        directed = m.captures[2] == "D", edges = Tuple{String,String}[]))
            continue
        end
        m = match(r"^\{\"([^\"]+)\",\s*\"([^\"]+)\"\},?$", line)
        m === nothing && error("unrecognized line in $(path): $(raw)")
        isempty(out) && error("edge before first system header in $(path)")
        push!(out[end].edges, (String(m.captures[1]), String(m.captures[2])))
    end
    return out
end
