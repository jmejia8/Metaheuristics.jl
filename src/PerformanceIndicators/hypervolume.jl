module HyperVolume
#=
Based on the file:
http://ls11-www.cs.uni-dortmund.de/_media/rudolph/hypervolume/hv_python.zip
written in python by Simon Wessing (2010).

This implementation is for the variant 3 detailed in:
C. M. Fonseca, L. Paquete, and M. Lopez-Ibanez. An improved dimension-sweep
algorithm for the hypervolume indicator. In IEEE Congress on Evolutionary
Computation, pages 1157-1163, Vancouver, Canada, July 2006.

=#

include("hypervolume-utils.jl")

"""
Sets up the list data structure needed for calculation.
"""
function preProcess(front, referencePoint)
    dimensions = length(referencePoint)
    sentinel = Node(dimensions)
    sentinel.next = fill(sentinel, dimensions)
    sentinel.prev = fill(sentinel, dimensions)
    nodes = [Node(dimensions, point) for point in front]

    for i in 1:dimensions
        # sortByDimension i
        sort!(nodes, by = node -> node.cargo[i])
        extend!(sentinel,nodes, i)
    end

    return sentinel

end



"""Recursive call to hypervolume calculation.

In contrast to the paper, the code assumes that the reference point
is [0, ..., 0]. This allows the avoidance of a few operations.
"""
function hvRecursive(sentinel, dimIndex, len, bounds)
    hvol = 0.0
    if len == 0
        return hvol
    elseif dimIndex == 1
        # special case: only one dimension
        # why using hypervolume at all?
        return -sentinel.next[1].cargo[1]

    elseif dimIndex == 2
        # special case: two dimensions, end recursion
        q = sentinel.next[2]
        h = q.cargo[1]
        p = q.next[2]

        while !(p === sentinel) #p is not sentinel
            pCargo = p.cargo
            hvol += h * (q.cargo[2] - pCargo[2])
            if pCargo[1] < h
                h = pCargo[1]
            end
            q = p
            p = q.next[2]
        end
        hvol += h * q.cargo[2]
        return hvol
    else
        p = sentinel
        q = p.prev[dimIndex]

        while !isempty(q.cargo) #q.cargo != None
            if q.ignore < dimIndex-1
                q.ignore = 0
            end
            q = q.prev[dimIndex]
        end

        q = p.prev[dimIndex]
        while len > 1 && (
                          q.cargo[dimIndex] > bounds[dimIndex] || q.prev[dimIndex].cargo[dimIndex] >= bounds[dimIndex])
            p = q
            remove!(p, dimIndex, bounds)
            q = p.prev[dimIndex]
            len -= 1
        end

        qArea = q.area
        qCargo = q.cargo
        qPrevDimIndex = q.prev[dimIndex]

        if len > 1
            hvol = qPrevDimIndex.volume[dimIndex] 
            hvol += qPrevDimIndex.area[dimIndex] * ( qCargo[dimIndex] - qPrevDimIndex.cargo[dimIndex])
        else
            qArea[1] = 1
            qArea[2:dimIndex] = [qArea[i] * -qCargo[i] for i in 1:dimIndex-1]
        end

        q.volume[dimIndex] = hvol
        if q.ignore >= dimIndex-1
            qArea[dimIndex] = qPrevDimIndex.area[dimIndex]
        else

            qArea[dimIndex] = hvRecursive(sentinel, dimIndex - 1, len, bounds)
            if qArea[dimIndex] <= qPrevDimIndex.area[dimIndex]
                q.ignore = dimIndex-1
            end
        end

        while !(p === sentinel) #p is not sentinel
            pCargoDimIndex = p.cargo[dimIndex]
            hvol += q.area[dimIndex] * (pCargoDimIndex - q.cargo[dimIndex])
            bounds[dimIndex] = pCargoDimIndex
            reinsert!(p, dimIndex, bounds)
            len += 1
            q = p
            p = p.next[dimIndex]
            q.volume[dimIndex] = hvol

            if q.ignore >= dimIndex-1
                q.area[dimIndex] = q.prev[dimIndex].area[dimIndex]
            else
                q.area[dimIndex] = hvRecursive(sentinel, dimIndex - 1, len, bounds)
                if q.area[dimIndex] <= q.prev[dimIndex].area[dimIndex]
                    q.ignore = dimIndex-1
                end
            end

        end
        hvol -= q.area[dimIndex] * q.cargo[dimIndex]
        return hvol


    end


end

function hv(front, referencePoint)
    dimensions = length(referencePoint)

    relevantPoints = [fx - referencePoint for fx in front]
    sentinel = preProcess(relevantPoints, referencePoint)

    bounds = fill(-Inf, dimensions)
    
    return hvRecursive(sentinel, dimensions, length(relevantPoints), bounds)
end

end




"""
	hypervolume(front, reference_point)
Computes the hypervolume indicator, i.e., volume between points in `front` and `reference_point`.

Note that each point in `front` must (weakly) dominates to `reference_point`. Also, `front`
is a non-dominated set.

If `front::State` and `reference_point::Vector`, then computes `hypervolume(front.population, reference_point)` after
ignoring solutions in `front` that do not dominate `reference_point`.
"""
function hypervolume(front::Array{Vector{T}}, reference_point::Vector; verbose=true) where T <: Real

	weaklyDominates(point, other) = begin
		for i in 1:length(point)
			if point[i] > other[i]
				return false
			end
		end
		return true
	end


	relevantPoints = Vector[]
	for point in front
		# only consider points that dominate the reference point
		if weaklyDominates(point, reference_point)
			push!(relevantPoints, point)
		end
	end

	if length(relevantPoints) != length(front)
		ign = length(front) - length(relevantPoints)
		rel = length(relevantPoints)
		verbose && @warn "Ignoring $ign points dominated by the reference point ($rel points are used)."
	end

	return HyperVolume.hv(relevantPoints, reference_point)
end

hypervolume(front::Vector{T}, reference_point::Vector{T}; verbose=true) where T <: AbstractMultiObjectiveSolution = hypervolume(fval.(front), fval(reference_point); verbose)
hypervolume(front::Vector{T}, reference_point::Vector; verbose=true) where T <: AbstractMultiObjectiveSolution = hypervolume(fval.(front), reference_point;verbose)

function hypervolume(front::Matrix, reference_point::Vector; verbose=true)
	front_ = [ front[i,:] for i in 1:size(front,1) ]
	hypervolume(front_, reference_point; verbose)
end
