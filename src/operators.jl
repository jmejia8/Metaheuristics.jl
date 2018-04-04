function Selection(Old::xf_indiv, New::xf_indiv, searchType=:minimize)
    if searchType == :minimize
        return New.f < Old.f
    end
    
    return New.f > Old.f
end

# Deb rules (selection)
function Selection(Old::xfgh_indiv, New::xfgh_indiv, searchType=:minimize)

    old_vio = countViolations(Old.g, Old.h)
    new_vio = countViolations(New.g, New.h)

    if new_vio < old_vio 
        return true
    elseif new_vio > old_vio 
        return false
    end

    if searchType == :minimize
         return New.f < Old.f
    end
    
    return New.f > Old.f
end

# Deb rules (selection)
function Selection(Old::xfg_indiv, New::xfg_indiv, searchType=:minimize)

    old_vio = countViolations(Old.g, [])
    new_vio = countViolations(New.g, [])

    if new_vio < old_vio 
        return true
    elseif new_vio > old_vio 
        return false
    end

    if searchType == :minimize
         return New.f < Old.f
    end
    
    return New.f > Old.f
end

function getBest(Population, searchType::Symbol = :minimize)
    best = Population[1]

    for i = 2:length(Population)
        if Selection(best, Population[i])
            best = Population[i]
        end
    end

    return best
end

function getBestInd(Population, searchType::Symbol = :minimize)
    j = 1

    for i = 2:length(Population)
        if Selection(Population[j], Population[i])
            j = i
        end
    end

    return j
end

function generateChild(individual::Type{xf_indiv}, x::Vector{Float64}, fResult::Float64)
    return individual(x, fResult)
end

function generateChild(individual::Type{xfg_indiv}, x::Vector{Float64}, fResult)
    f, g = fResult
    return individual(x, f, g)
end

function generateChild(individual::Type{xfgh_indiv}, x::Vector{Float64}, fResult)
    f, g, h = fResult
    return individual(x, f, g, h)
end