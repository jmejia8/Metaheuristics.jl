#=
Based on the file:
http://ls11-www.cs.uni-dortmund.de/_media/rudolph/hypervolume/hv_python.zip
originally written in python by Simon Wessing (2010).

This implementation is for the variant 3 detailed in:
C. M. Fonseca, L. Paquete, and M. Lopez-Ibanez. An improved dimension-sweep
algorithm for the hypervolume indicator. In IEEE Congress on Evolutionary
Computation, pages 1157-1163, Vancouver, Canada, July 2006.
=#

mutable struct Node
    cargo::Vector
    next::Vector{Union{Nothing,Node}}
    prev::Vector{Union{Nothing,Node}}
    ignore::Int
    area::Vector{Float64}
    volume::Vector{Float64}
end

function Node(numberLists, cargo = zeros(0);
        next   = fill(nothing, numberLists),
        prev   = fill(nothing, numberLists),
        ignore = 0,
        area   = zeros(numberLists),
        volume = zeros(numberLists),
    )

    Node(cargo, next, prev, ignore, area, volume)
    
end

function extend!(sentinel::Node, nodes, index)
    for node in nodes
        lastButOne = sentinel.prev[index]
        node.next[index] = sentinel
        node.prev[index] = lastButOne
        # set the last element as the new one
        sentinel.prev[index] = node
        lastButOne.next[index] = node
    end
end



"""
    Removes and returns 'node' from all lists in 1,2,..., 'index'.
"""
function remove!(node::Node, index, bounds)
    for i in 1:index-1
        predecessor = node.prev[i]
        successor = node.next[i]
        predecessor.next[i] = successor
        successor.prev[i] = predecessor

        if bounds[i] > node.cargo[i]
            bounds[i] = node.cargo[i]
        end
    end
    return node
end


"""
Inserts 'node' at the position it had in all lists in [1, 'index']
before it was removed. This method assumes that the next and previous
nodes of the node that is reinserted are in the list.

"""
function reinsert!(node::Node, index, bounds)
    for i in 1:index
        node.prev[i].next[i] = node
        node.next[i].prev[i] = node
        if bounds[i] > node.cargo[i]
            bounds[i] = node.cargo[i]
        end
    end

    node
end


