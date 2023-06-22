mutable struct InsertRandomMutation
    num_mutants::Int
end

function mutation!(Q, m::InsertRandomMutation)
    # TODO: add rng
    R = rand(m.num_mutants, size(Q, 2))
    vcat(Q, R)
end
