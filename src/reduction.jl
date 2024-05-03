import AbstractAlgebra as AA

function pseudo_reduce(f::RE, G) where {RE<:AA.MPolyRingElem}
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
            f_ = AA.leading_coefficient(g)*f_ - a*γ*g
        end
    end
    return c, r
end

function inter_reduce(G::Vector{RE}) where {RE<:AA.MPolyRingElem}
    for i in 1:length(G)
        c, r = pseudo_reduce(G[i], G[[1:(i-1) ; (i+1):end]])
        G[i] = r
    end
    filter(g -> !iszero(g), G)
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
