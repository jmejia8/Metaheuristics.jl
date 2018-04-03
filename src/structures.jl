struct xf_indiv # Single Objective
	x::Vector{Float64}
	f::Float64
end

struct xfg_indiv # Single Objective Constraied
	x::Vector{Float64}
	f::Float64
	g::Vector{Float64}
end

struct xfgh_indiv # Single Objective Constraied
	x::Vector{Float64}
	f::Float64
	g::Vector{Float64}
	h::Vector{Float64}
end

struct xFgh_indiv # Single Objective Constraied
	x::Vector{Float64}
	f::Vector{Float64}
	g::Vector{Float64}
	h::Vector{Float64}
end