import JMcDM: mcdm, MCDMSetting, summary, makeminmax,DataFrame, MCDMMethod

export mcdm
export MCDMSetting
export summary

export best_alternative

function JMcDM.mcdm(
    f::AbstractMatrix{<:AbstractFloat}, # objective functions by col
    w::AbstractVector{<:AbstractFloat},
    method::T
    ) where T <: MCDMMethod

    fns = makeminmax([minimum for i in 1:size(f, 1)])
    mcdm(DataFrame(f, :auto), w, fns, method)
end

function JMcDM.mcdm(population::AbstractArray{<: AbstractMultiObjectiveSolution}, args...)
    fs = fvals(population) # Pareto solutions
    mcdm(fs, args...)
end

function JMcDM.mcdm(st::State,args...)
    mcdm(st.population, args...)
end

function JMcDM.MCDMSetting(f::AbstractMatrix{<: AbstractFloat}, weights)
    fns = makeminmax([minimum for i in 1:size(f, 1)])
    MCDMSetting(DataFrame(f, :auto), weights, fns)
end


function JMcDM.MCDMSetting(status::State,args...)
    MCDMSetting(status.population, args...)
end

function JMcDM.MCDMSetting(population::AbstractArray{<: AbstractMultiObjectiveSolution}, args...)
    MCDMSetting(fvals(population), args...)
end

function JMcDM.summary(
    f::AbstractMatrix{<: AbstractFloat},
    w::Array{Float64,1}, 
    methods::Array{Symbol,1}
    )

    fns = makeminmax([minimum for i in 1:size(f, 1)]) 
    summary(DataFrame(f, :auto), w, fns, methods)
end


function JMcDM.summary(population::AbstractArray{<: AbstractMultiObjectiveSolution},args...)
    summary(fvals(population), args...)
end


function JMcDM.summary(st::State,args...)
    summary(st.population, args...)
end

# const decisionmaking = mcdm

function best_alternative(
    population::AbstractArray{<: AbstractMultiObjectiveSolution},
    w::AbstractVector{<:AbstractFloat},
    method::T
    ) where T <: MCDMMethod
    result = mcdm(population, w, method)
    population[result.bestIndex]
end


function best_alternative(population::AbstractArray{<: AbstractMultiObjectiveSolution}, args...)
    best_alternative(fvals(population), args...)
end


function best_alternative(st::State, args...)
    best_alternative(st.population, args...)
end

best_alternative(status::State, args...) = best_alternative(status.population,args...)

