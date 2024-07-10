import AbstractAlgebra as AA
import AbstractAlgebra.Generic.MPoly

@doc raw"""
    pseudo_reduce(f::RE, G, reduced = false) where {RE<:MPoly}

Pseudo-reduces ``f`` modulo ``G = \{g_1, g_2, \dots, g_n\}``. Pseudo-reduction writes
```math
c f = r + \sum_{i=1}^n f_i g_i
```
where ``c`` is a product of leading coefficients of ``G``,
``\operatorname{lm}(f_i g_i) \leq \operatorname{lm}(f)``,
and no term of ``r`` is divisible by any leading monomial of ``G``.

The result of `pseudo_reduce(f, G)` is the pair `(c, r)`.

Takes the argument `reduced`. If set to `true`, the function will remove factors in `r`
 coming from  leading coefficients of ``G``.
"""
function pseudo_reduce(f::RE, G, reduced = false) where {RE<:MPoly}
    AX = AA.parent(f)
    A = AA.base_ring(AX)
    r = zero(AX)
    c = one(A)
    f_ = f

    while !iszero(f_)
        a, m = AA.leading_coefficient(f_), AA.leading_monomial(f_)
        divs = [AA.divides(m, AA.leading_monomial(g)) for g in G]
        i = findfirst(d -> d[1], divs)
        if isnothing(i)
            r = r + a*m
            f_ = AA.tail(f_)
        else
            g = G[i]
            γ = divs[i][2]
            c = c*AA.leading_coefficient(g)
            r = r*AA.leading_coefficient(g)
            f_ = AA.leading_coefficient(g)*f_ - a*γ*g
        end
    end

    if reduced
        H = AA.leading_coefficient.(G)
        for h ∈ H
            i = 0
            while AA.divides(r, (h^(i+1))*one(AX))[1]
                i = i + 1
            end
            r = r / h^i
        end
    end
    return c, r
end

"""
    pseudo_remainder(f::RE, G, reduced = false) where {RE<:MPoly}

Defined as `pseudo_reduce(f, G, reduced)[2]`.

See also [`pseudo_reduce`](@ref)
"""
function pseudo_remainder(f::RE, G, reduced = false) where {RE<:MPoly}
    return pseudo_reduce(f, G, reduced)[2]
end


function faithful_pseudo_reduce(f::Tuple{RE, RE}, G) where {RE<:MPoly}
    AX = AA.parent(f[1])
    A = AA.base_ring(AX)
    r = zero(AX)
    c = one(A)
    f_ = f

    while !iszero(f_[1])
        a, m = AA.leading_coefficient(f_[1]), AA.leading_monomial(f_[1])
        divs = [AA.divides(m, AA.leading_monomial(g[1])) for g in G]
        i = findfirst(d -> d[1], divs)
        if isnothing(i)
            r = r + a*m
            f_ = f_ .- a*m
        else
            g = G[i]
            γ = divs[i][2]
            c = c*AA.leading_coefficient(g[1])
            r = r*c
            f_ = AA.leading_coefficient(g[1]).*f_ .- a*γ.*g
        end
    end
    return c, r, r + f_[2]
end


function inter_reduce(G::Vector{RE}) where {RE<:MPoly}
    G_ = empty(G)
    for g in G
        if !any([AA.divides(AA.leading_monomial(g), AA.leading_monomial(g_))[1] for g_ in G_])
            push!(G_, g)
        end
    end
    for i in 1:length(G_)
        c, r = pseudo_reduce(G_[i], G_[[1:(i-1) ; (i+1):end]])
        G_[i] = r
    end
    filter(g -> !iszero(g), G_)
end

function faithful_inter_reduce(G::Vector{Tuple{RE, T}}) where {RE<:MPoly, T}
    A = AA.base_ring(G[1][1])
    G_ = empty(G)
    for (g, g_f) in G
        if !any([AA.divides(AA.leading_monomial(g), AA.leading_monomial(g_[1]))[1] for g_ in G_])
            push!(G_, (g, g_f))
        end
    end
    G__ = empty([(G[1][1], G[1][2])])
    for i in 1:length(G_)
        c, r, r_f = faithful_pseudo_reduce(G_[i], [g_ for g_ in G_[[1:(i-1) ; (i+1):end]]])
        push!(G__, (r, r_f))
    end
    filter(g -> !iszero(g[2]), G__)
end




function pseudo_reduce(f, G, params)
    f_ = make_parameters(f, params)
    G_ = [make_parameters(g, params) for g in G]
    r = zero(f)
    c = one(f)

    while !iszero(f_)
        t = leading_term(f_)
        a, m = coefficient(t), monomial(t)
        i = findfirst(g -> divides(leading_monomial(g), m), G_)
        if isnothing(i)
            r = r + t
            f_ = f_ - t
        else
            g = G_[i]
            c = c*leading_coefficient(g)
            γ = div(m, leading_monomial(g))
            f_ = term(leading_coefficient(g), one(m))*f_ - term(a, one(m))*γ*g
        end
    end
    return c, unparameterize(r)
end

function inter_reduce(G, params)
    for i in 1:length(G)
        c, r = pseudo_reduce(G[i], G[[1:(i-1) ; (i+1):end]], params)
        G[i] = r
    end
    filter(g -> !iszero(g), G)
end



export pseudo_reduce, inter_reduce
