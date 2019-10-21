#
# Test routines for the module SparseMaxFlowMinCut
#
# by Artur Pessoa (2019)
#

include("../src/SparseMaxFlowMinCut.jl")

function run_test(inst::Int)
    # check the number of program arguments
    if (length(ARGS) > 1)
        println("Use: julia runtests.jl [<filename>]")
        exit(-1)
    end

    # read all the input file to the vector "numbers"
    fname = "data/inst$(inst).txt"
    if !isempty(ARGS)
        fname = ARGS[1]
    end
    f = open(fname, "r")
    numbers = map(y->parse(Int,y),split(read(f,String)))
    close(f)

    # get the data from "numbers"
    n = Int(numbers[1])
    m = Int(numbers[2])
    s = Int(numbers[3])
    t = Int(numbers[4])
    pos = 5
    A = SparseMaxFlowMinCut.ArcFlow[]
    for a = 1:m
        push!(A, SparseMaxFlowMinCut.ArcFlow(numbers[pos], numbers[pos + 1], numbers[pos + 2]))
        pos += 3
    end
    check = numbers[pos]
    if isempty(ARGS)
        println("inst$(inst): $check")
    else
        @show n, m, s, t, A
    end

    # call the maximum flow/minimum cut function
    maxFlow, flows, cut = SparseMaxFlowMinCut.find_maxflow_mincut(
        SparseMaxFlowMinCut.Graph(n, A), s, t
    )
    @assert check == maxFlow
    return maxFlow, flows, cut
end

if !isempty(ARGS)
    @time res = run_test(0)
    @show res
else
    for i in 0:15
        run_test(i)
    end
    println("Ok!")
end
