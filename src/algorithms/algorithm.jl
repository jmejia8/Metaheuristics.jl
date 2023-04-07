iscompatible(search_space, algorithm) = false
iscompatible(search_space::BoxConstrainedSpace, algorithm::AbstractParameters) = true
iscompatible(search_space::AbstractMatrix, algorithm::AbstractParameters) = true

function check_compat(search_space, algorithm::AbstractParameters)
    if iscompatible(search_space, algorithm)
        return
    end
    
    a = typeof(algorithm)
    s = typeof(search_space)
    error("Implemented metaheuristic `$a` is not compatible with search space `$s`.")
end

