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
