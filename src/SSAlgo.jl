using Groebner

function is_parameters(g::AP, params) where {AP<:AbstractPolynomial}
    return maxdegree(make_parameters(g, params)) == 0
end

function CGSMain(F::Vector{AP}, params) where {AP<:AbstractPolynomial}
    G = groebner(F)
    if one(F[1]) ∈ G
        return [(one(F[1]), G)] # Skal det være F eller G til sidst?
    else
        println(G)
        hs = [leading_coefficient(make_parameters(g, params)) for g ∈ G
              if !is_parameters(g, params)]
        h = reduce(lcm, hs)
        if h == one(h)
            return [(h, G)]
        else
            return vcat([(h, G)], [CGSMain([G ; [hi]], params) for hi in hs]...)
        end
    end
end

function CGS(F, params)
    H = CGSMain(F, params)
    G₀ = groebner(F)
    GS = empty([(F, F, F)])
    if any([is_parameters(g, params) for g ∈ G₀])
        GS = [(empty(F), filter(g -> is_parameters(g, params), G₀), [one(F[1])])]
    end
    for (h, G) ∈ H
        push!(GS, (filter(g -> is_parameters(g, params), G), [h], filter(g -> !is_parameters(g, params), G)))
    end
    return GS
end

function CGBMain(F, S, params)
    if one(F[1]) ∈ groebner(S)
        return empty([(F, F, F)])
    else
        vars = collect(Set(vcat([variables(p) for p ∈ F]...)))
        @polyvar __x[1:length(vars) + 1] # Bem. __x[1] > __x[2] > ... Det holder både i Lex og RevLex
        F = [subs(f, (v => x for (v, x) ∈ zip(vars, __x[2:end]))...) for f ∈ F]
        S = [subs(s, (v => x for (v, x) ∈ zip(vars, __x[2:end]))...) for f ∈ F]

        G = groebner([[__x[1]*f for f ∈ F] ; [(__x[1]-1)*s for s ∈ S]])
    end
end
