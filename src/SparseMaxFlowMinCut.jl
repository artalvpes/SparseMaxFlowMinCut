module SparseMaxFlowMinCut

using DataStructures

#
# Julia implementation of a shortest augmenting path algorithm for the
# maximum flow/minimu cut problem using sparse graphs
#
# by Artur Pessoa (2019)
#

export ArcFlow, Graph, find_maxflow_mincut

struct ArcFlow
    i::Int
    j::Int
    f::Int
end

mutable struct Graph
    _n::Int
    adj::Vector{Vector{Int}}
    val::Vector{Vector{Int}}

    # position of the opposite arc in the adjacency list of its head node
    opp::Vector{Vector{Int}}
end

# constructor
function Graph(n::Int, A::Vector{ArcFlow})
    g = Graph(n, [Vector{Int}() for i=1:n], [Vector{Int}() for i=1:n],
        [Vector{Int}() for i=1:n]
    )

    # add the arcs to the graph
    for a in A
        push!(g.adj[a.i], a.j)
        push!(g.val[a.i], a.f)
        push!(g.opp[a.i], length(g.opp[a.j])+1)
        push!(g.adj[a.j], a.i)
        push!(g.val[a.j], 0)
        push!(g.opp[a.j], length(g.opp[a.i]))
    end
    
    return g
end

# reset the arc values in the graph
function reset(g::Graph, A::Vector{ArcFlow}, buffer::Vector{Int})
    arc_pos = buffer
    resize!(arc_pos, g._n)
    fill!(arc_pos, 0)
    for a in A
        arc_pos[a.i] += 1
        g.val[a.i][arc_pos[a.i]] = a.f
        arc_pos[a.j] += 1
        g.val[a.j][arc_pos[a.j]] = 0
    end
end

# find a path p from s to t using only non-zero arcs and return its minimum arc value, and p
# * the path is represented by a sequence of vertex indices in adjacency lists starting at
#   t until reach s.
# * if no such a path exist, store in p a list of vertices that are reachable from s and
#   return 0, and p.
function find_path(g::Graph, s::Int, t::Int)::Tuple{Int, Vector{Int}}
    # do a bfs in the graph
    prev = fill(-2, g._n)
    q = CircularDeque{Int}(g._n)
    push!(q, s)
    prev[s] = -1
    while !isempty(q)
        i = popfirst!(q)

        # check if the sink was reached
        if i == t
            # fill the path and return the minimum value
            minVal = typemax(Int)
            p = Vector{Int}()
            while prev[i] != -1
                push!(p, prev[i])
                k = g.opp[i][prev[i]]
                i = g.adj[i][prev[i]]
                if (minVal > g.val[i][k])
                    minVal = g.val[i][k]
                end
            end
            return minVal, p
        end

        # check the adjacent vertices
        l = length(g.adj[i])
        for k in 1:l
            if g.val[i][k] == 0
                continue
            end
            j = g.adj[i][k]
            if prev[j] != -2
                continue
            end
            prev[j] = g.opp[i][k]
            push!(q, j)
        end
    end

    # store in p a list of vertices that are reachable from s and return 0
    p = Vector{Int}()
    for i = 1:g._n
        if prev[i] != -2
            push!(p, i)
        end
    end
    return 0, p
end

# find a minimum s-t cut in the graph instance and return the cut value
function find_maxflow_mincut(instance::Graph, s::Int, t::Int
        )::Tuple{Int, Vector{ArcFlow}, Vector{Bool}}
    if s == t
        return typemax(Int), ArcFlow[], Bool[]
    end

    # make a copy of the instance graph as the residual graph
    residual = deepcopy(instance)

    # iterate finding paths from s to t
    p = Vector{Int}()
    flow, p = find_path(residual, s, t)
    maxFlow = 0
    while (flow > 0)
        # pass flow through the path
        i = t
        ll = length(p)
        for l in 1:ll
            ki = p[l]
            j = residual.adj[i][ki]
            kj = residual.opp[i][ki]
            residual.val[j][kj] -= flow
            residual.val[i][ki] += flow
            i = j
        end
        @assert i == s
        maxFlow += flow
        flow, p = find_path(residual, s, t)
    end
    
    # build the vector of arc flows
    flows = ArcFlow[]
    for i in 1:instance._n
        l = length(instance.adj[i])
        for ll in 1:l
            if instance.val[i][ll] > 0
                push!(flows, ArcFlow(i, instance.adj[i][ll],
                    instance.val[i][ll] - residual.val[i][ll])
                )
            end
        end
    end

    # here, p contains a list of vertices in the cut... convert it the boolean vector cut
    cut = zeros(Bool, residual._n)
    ll = length(p)
    for l in 1:ll
        cut[p[l]] = true
    end
    return maxFlow, flows, cut
end

end # module
