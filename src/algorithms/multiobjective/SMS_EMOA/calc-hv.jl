# Exclusive hypervolume contribution of each point (column) of `P` (M × N matrix) inside
# the box bounded by `ref`. Points are assumed mutually non-dominated and strictly below
# `ref`.

"""
    hv_contrib_2d(P)

Exact contributions for two objectives. Extreme points get `Inf`.
"""
function hv_contrib_2d(P::AbstractMatrix{<:Real})
    N = size(P, 2)
    ΔS = fill(Inf, N)
    N < 3 && return ΔS
    # in a non-dominated front, sorting by f₁ sorts f₂ in reverse order
    p = sortperm(1:N, by = i -> (P[1, i], P[2, i]))
    @inbounds for t in 2:N - 1
        i = p[t]
        ΔS[i] = (P[1, p[t + 1]] - P[1, i]) * (P[2, p[t - 1]] - P[2, i])
    end
    return ΔS
end

"""
    hv_contrib_3d(P, ref)

Exact contributions for three objectives. The contribution of `p` is the volume of the box
`[p, ref]` minus the hypervolume of `{max(q, p) : q ≠ p}` inside that box.
"""
function hv_contrib_3d(P::AbstractMatrix{<:Real}, ref)
    N = size(P, 2)
    T = Matrix{Float64}(undef, 3, N) # other points projected on the box of p
    xs = Vector{Float64}(undef, N)   # 2D staircase used by the sweep
    ys = Vector{Float64}(undef, N)
    perm = Vector{Int}(undef, N)
    r1, r2, r3 = ref[1], ref[2], ref[3]
    ΔS = zeros(N)

    @inbounds for i in 1:N
        a, b, c = P[1, i], P[2, i], P[3, i]
        n = 0
        for q in 1:N
            q == i && continue
            x = max(P[1, q], a)
            y = max(P[2, q], b)
            z = max(P[3, q], c)
            (x < r1 && y < r2 && z < r3) || continue
            n += 1
            T[1, n] = x
            T[2, n] = y
            T[3, n] = z
        end
        box = (r1 - a) * (r2 - b) * (r3 - c)
        ΔS[i] = max(0.0, box - _hv3d!(T, n, ref, xs, ys, perm))
    end
    return ΔS
end

# Hypervolume of the first `n` columns of `T` (3 × n) w.r.t. `ref`: sweep on f₃ keeping the
# 2D non-dominated staircase of (f₁, f₂) (xs ascending, ys descending) and its area.
function _hv3d!(T, n, ref, xs, ys, perm)
    n == 0 && return 0.0
    r1, r2, r3 = ref[1], ref[2], ref[3]
    p = view(perm, 1:n)
    sortperm!(p, view(T, 3, 1:n), alg = Base.Sort.QuickSort) # no scratch allocation
    k = 0
    area = 0.0
    vol = 0.0
    zprev = T[3, p[1]]
    @inbounds for t in 1:n
        i = p[t]
        x, y, z = T[1, i], T[2, i], T[3, i]
        vol += area * (z - zprev)
        zprev = z

        # first position with xs ≥ x
        j = searchsortedfirst(view(xs, 1:k), x)
        # dominated by the previous step or by an equal x
        (j > 1 && ys[j - 1] <= y) && continue
        (j <= k && xs[j] == x && ys[j] <= y) && continue

        # points dominated by (x, y) are contiguous from j
        e = j
        while e <= k && ys[e] >= y
            e += 1
        end
        xe = e <= k ? xs[e] : r1

        # area added by (x, y)
        h = j > 1 ? ys[j - 1] : r2
        area += ((j <= k ? xs[j] : r1) - x) * (h - y)
        for m in j:e - 1
            area += ((m + 1 < e ? xs[m + 1] : xe) - xs[m]) * (ys[m] - y)
        end

        # replace xs[j:e-1] by (x, y)
        nrem = e - j
        if nrem == 0
            for m in k:-1:j
                xs[m + 1] = xs[m]
                ys[m + 1] = ys[m]
            end
        elseif nrem > 1
            for m in e:k
                xs[m - nrem + 1] = xs[m]
                ys[m - nrem + 1] = ys[m]
            end
        end
        k += 1 - nrem
        xs[j] = x
        ys[j] = y
    end
    vol += area * (r3 - zprev)
    return vol
end

"""
    hv_contrib_mc(P, ref, n_samples, rng)

Monte Carlo estimate of the exclusive contributions (HypE with k = 1), for any number of
objectives. Samples are drawn uniformly in `[min(P), ref]`; a sample counts for a point
only when that point is its only dominator. One pass over samples, O(N) memory.
"""
function hv_contrib_mc(P::AbstractMatrix{<:Real}, ref, n_samples, rng = default_rng_mh())
    M, N = size(P)
    lo = vec(minimum(P, dims = 2))
    w = ref .- lo
    hits = zeros(Int, N)
    s = Vector{Float64}(undef, M)
    @inbounds for _ in 1:n_samples
        for m in 1:M
            s[m] = lo[m] + w[m] * rand(rng)
        end
        owner = 0
        for i in 1:N
            dom = true
            for m in 1:M
                if P[m, i] > s[m]
                    dom = false
                    break
                end
            end
            dom || continue
            if owner == 0
                owner = i
            else
                # dominated by two or more points
                owner = -1
                break
            end
        end
        owner > 0 && (hits[owner] += 1)
    end
    return hits .* (prod(w) / n_samples)
end
