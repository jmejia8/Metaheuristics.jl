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

    old_vio = Old.sum_violations
    new_vio = New.sum_violations

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
    old_vio = Old.sum_violations
    new_vio = New.sum_violations

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


