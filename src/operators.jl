function Selection(Old::xf_indiv, New::xf_indiv, searchType::Symbol=:minimize; leq::Bool=false)
    if searchType == :minimize
        if leq
            return New.f <= Old.f
        end

        return New.f < Old.f
    end

    if leq
        return New.f >= Old.f
    end

    return New.f > Old.f
end

# Deb rules (selection)
function Selection(Old::xfgh_indiv, New::xfgh_indiv, searchType::Symbol=:minimize; leq::Bool=false)

    old_vio = violationsSum(Old.g, Old.h)
    new_vio = violationsSum(New.g, New.h)

    if new_vio < old_vio
        return true
    elseif new_vio > old_vio
        return false
    end

    if searchType == :minimize
        if leq
            return New.f <= Old.f
        end
        return New.f < Old.f
    end

    if leq
        return New.f >= Old.f
    end

    return New.f > Old.f
end

# Deb rules (selection)
function Selection(Old::xfg_indiv, New::xfg_indiv, searchType::Symbol=:minimize; leq::Bool=false)
    old_vio = violationsSum(Old.g, [])
    new_vio = violationsSum(New.g, [])

    if new_vio < old_vio
        return true
    elseif new_vio > old_vio
        return false
    end

    if searchType == :minimize
        if leq
            return New.f <= Old.f
        end
        return New.f < Old.f
    end
    if leq
        return New.f >= Old.f
    end

    return New.f > Old.f
end

function getBest(Population, searchType::Symbol = :minimize, is_better=is_better)
    best = Population[1]

    for i = 2:length(Population)
        if is_better(Population[i], best)
            best = Population[i]
        end
    end

    return best
end

function getWorst(Population, searchType::Symbol = :minimize, is_better=is_better)
    worst = Population[1]

    for i = 2:length(Population)
        if is_better(Population[i], worst)
            worst = Population[i]
        end
    end

    return worst
end

function getWorstInd(Population, searchType::Symbol = :minimize, is_better=is_better)
    worst = 1

    for i = 2:length(Population)
        if is_better(Population[worst], Population[i])
            worst = i
        end
    end

    return worst
end

function is_better(x, y)
    # x better than y
    return Selection(y, x)
end

function getBestInd(Population, searchType::Symbol = :minimize, is_better = is_better)
    j = 1

    for i = 2:length(Population)
        if is_better(Population[i], Population[j])
            j = i
        end
    end

    return j
end

function generateChild(x::Vector{Float64}, fResult::Float64)
    return xf_indiv(x, fResult)
end

function generateChild(x::Vector{Float64}, fResult::Tuple{Float64,Array{Float64,1}})
    f, g = fResult
    return xfgh_indiv(x, f, g, [0.0])
end

function generateChild(x::Vector{Float64}, fResult::Tuple{Float64,Array{Float64,1},Array{Float64,1}})
    f, g, h = fResult
    return xfgh_indiv(x, f, g, h)
end

function generateChild(x::Vector{Float64}, fResult::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}})
    f, g, h = fResult
    return xFgh_indiv(x, f, g, h)
end

function inferType(fVal::Tuple{Float64})
    return xf_indiv
end

function inferType(fVal::Tuple{Float64,Array{Float64,1}})
    return xfg_indiv
end

function inferType(fVal::Tuple{Float64,Array{Float64,1},Array{Float64,1}})
    return xfgh_indiv
end

function inferType(fVal::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}})
    return xFgh_indiv
end
