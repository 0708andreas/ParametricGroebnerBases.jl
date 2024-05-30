using Groebner
import MultivariatePolynomials.leading_coefficient,
    MultivariatePolynomials.ordering,
    MultivariatePolynomials.leading_monomial
import Base.lcm
# using AbstractAlgebra

function ordering(f::DynamicPolynomials.Polynomial{C, O, T}) where {C, O, T}
    return O
end

# function has_only(g::AP, params) where {AP<:AbstractPolynomial}
#     return isconstant(make_parameters(g, params))
# end

function has_only(g::AP, params) where {AP<:AbstractPolynomial}
    return all(effective_variables(g) .∈ Ref(params))
end

export has_only

function leading_coefficient(g::AP, params::Tuple) where {AP<:AbstractPolynomial}
    return leading_coefficient(make_parameters(g, params))
end

function leading_coefficient2(g::AP, params::Tuple) where {AP<:AbstractPolynomial}
    res = zero(g)
    s = (p => 1 for p ∈ variables(g) if p ∈ params)
    lm = monomial(subs(leading_term(g), s...))
    for t in reverse(terms(g))
        if lm == monomial(subs(t, s...))
            res = res + t
        else
            break
        end
    end
    return subs(res, (x => 1 for x ∈ variables(g) if !(x ∈ params))...)
end


function leading_monomial(g::AP, params::Tuple) where {AP<:AbstractPolynomial}
    return MultivariatePolynomials.leading_monomial(make_parameters(g, params))
end

function lcm(f::AP, g::AP) where {AP<:AbstractPolynomial{Int64}}
    return map_coefficients(c -> Int64(c), f * div(g, gcd(f, g)))
end

import AbstractAlgebra as AA

function parametrize(f, vars, params, AX)
    f_ = zero(AX)
    for (c, v) in zip(AA.coefficients(f), AA.exponent_vectors(f))
        f_ = f_ + (c * prod(vars .^ v[1:length(vars)]) * prod(params .^ v[(length(vars)+1):end]))
    end
    f_
end


function groebner_wrapper(F::Vector{RE}) where {RE <: AA.MPolyRingElem}
    if isempty(F) return F end

    AX = AA.parent(F[1])
    A = AA.base_ring(AX)
    K = AA.base_ring(A)
    params = AA.gens(A)
    params_symb = AA.symbols(A)
    vars = AA.gens(AX)
    vars_symb = AA.symbols(AX)

    KXY, news = AA.polynomial_ring(K, [vars_symb ; params_symb])
    newvars = news[1:length(vars)]
    newparams = news[length(vars)+1:end]

    vars_order = Dict(:lex => Lex(newvars),
                      :deglex => DegLex(newvars),
                      :degrevlex => DegRevLex(newvars))[AA.internal_ordering(AX)]
    order = ProductOrdering(vars_order, DegRevLex(newparams))

    F_ = [AA.evaluate(AA.map_coefficients(c -> AA.evaluate(c, newparams), f), newvars) for f ∈ F]
    G_ = groebner(F_, ordering=order)

    [parametrize(g, vars, params, AX) for g ∈ G_]
end

export groebner_wrapper

function CGS_simple(F::Vector{RE}) where {RE <: AA.MPolyRingElem}
    AX = AA.parent(F[1])
    A = AA.base_ring(AX)
    function CGS_inner(S)
        if !isempty(S)
            S = groebner(S)
        end
        if 1 ∈ S
            return empty([(F, F, F)])
        end

        G = groebner_wrapper([F ; one(AX) .* S])
        G = [g for g ∈ G if !iszero(divrem(AA.leading_coefficient(g), S)[2])]

        hs = unique([AA.leading_coefficient(g) for g ∈ G])
        h = one(AX)*reduce(lcm, hs, init=one(A))

        return vcat([(one(AX).*S, [h], G)],
                    [CGS_inner([S ; hi]) for hi ∈ hs]...)
    end
    return unique(CGS_inner(empty([one(A)])))
end


function CGS(F::Vector{RE}, red=true) where {RE <: AA.MPolyRingElem}
    AX = AA.parent(F[1])
    A = AA.base_ring(AX)
    function CGS_inner(S)
        if !isempty(S)
            S = groebner(S)
        end
        if 1 ∈ S
            return empty([(F, F, F)])
        end

        G = groebner_wrapper([F ; one(AX).*S])
        G = [g for g ∈ G if !iszero(divrem(AA.leading_coefficient(g), S)[2])]
        if red && !isempty(G)
            G = [G[i] for i in 1:length(G) if
                   !any([AA.divides(AA.leading_monomial(G[i]), AA.leading_monomial(g_))[1]
                         for g_ in G[1:i-1]]) ]
        end

        hs = unique([AA.leading_coefficient(g) for g ∈ G ])
        h = one(AX)*reduce(lcm, hs, init=one(A))

        if isempty(G)
            G = [zero(F[1])]
        end
        if red
            G = inter_reduce(filter(g -> !iszero(g), G))
        end

        return vcat([(one(AX).*S, [h], G)],
                    [CGS_inner([S ; [hi]]) for hi ∈ hs]...)
    end

    S = [AA.leading_coefficient(g) for g ∈ groebner_wrapper(F) if AA.is_constant(g)]
    if isempty(S)
        return unique(CGS_inner(S))
    else
        return unique([[(empty(F), one(AX).*S, [one(F[1])])] ; CGS_inner(S)])
    end
end


function CGS(F::Vector{AP}, params) where {AP<:AbstractPolynomial}
    vars = [v for v in variables(F) if v ∉ params]
    pars = collect(params)
    order = ProductOrdering(DegLex(vars), DegLex(pars))
    function CGS_inner(S)
        if 1 ∈ S || (!isempty(S) && 1 ∈ groebner(S))
            return empty([(F, F, F)])
        end

        G = groebner([F ; S], ordering=order)

        # @assert isgroebner(G, ordering=Lex())
        hs = unique([leading_coefficient(g, params) for g ∈ G
                  if !has_only(g, params)])
        h = reduce(lcm, hs, init=one(F[1]))

        if isempty(S)
            S_ = S
        else
            S_ = groebner(S)
        end

        G_ = [g for g ∈ G if !has_only(g, params) || isconstant(g)]
        if isempty(G_)
            G_ = [zero(F[1])]
        end

        return vcat([(S_, [h], G_)],
                    [CGS_inner([S ; [hi]]) for hi ∈ hs]...)
        # return vcat([(S_, [h], G)], [CGS_inner([S ; [hi]]) for hi ∈ hs]...)
    end

    S = [g for g ∈ groebner(F, ordering=order) if has_only(g, params)]
    if isempty(S)
        return CGS_inner(S)
    else
        return [[(empty(F), S, [one(F[1])])] ; CGS_inner(S)]
    end
end


export CGS

function CGBMain(F::Vector{RE1}, S::Vector{RE2}, red = true) where {RE1<:AA.MPolyRingElem, RE2<:AA.MPolyRingElem}
    if !isempty(S)
        S = groebner(S)
    end
    if 1 ∈ S
        return empty([(F, F, F)])
    end
    AX = AA.parent(F[1])
    vars = AA.gens(AX)
    vars_symb = AA.symbols(AX)
    A = AA.base_ring(AX)
    params = AA.gens(A)
    AtX, tvars = A[:t, vars_symb...]
    t = tvars[1]

    σ₁(g) = AA.evaluate(g, [[one(AX)] ; vars])
    σ₀(g) = AA.evaluate(g, [[zero(AX)] ; vars])

    F_ = [AA.evaluate(f, tvars[2:end]) for f ∈ F]
    G = groebner_wrapper([(t.*F_) ; ((1 - t).*S)])
    G = [g for g ∈ G if !iszero(divrem(AA.leading_coefficient(g), S)[2]) && AA.divides(AA.leading_term(g), t)[1]]
    if red && !isempty(G)
        G = [G[i] for i in 1:length(G)
                 if !any([AA.divides(AA.leading_monomial(G[i]), AA.leading_monomial(g_))[1]
                          for g_ in G[1:i-1]])]
    end

    hs = unique([AA.leading_coefficient(g) for g ∈ G
                     if AA.divides(AA.leading_term(g), t)[1] ])
    h = one(AX)*reduce(lcm, hs, init=one(A))

    G_ = [(g - g_t, g) for (g, g_t) ∈ zip(σ₁.(G), σ₀.(G)) ]
    if red && !isempty(G_)
        G_ = faithful_inter_reduce(G_)
        G_ = [g_f for (g_, g_f) in G_]
    else
        G_ = [g_f for (g, g_f) ∈ G_]
    end

    return vcat([(one(AX).*S, [h], G_)],
                [CGBMain(F, [S ; [hi]], red) for hi in hs]...)
end

function CGB(F::Vector{RE}, red=true) where {RE<:AA.MPolyRingElem}
    # S = [AA.leading_coefficient(g) for g ∈ groebner_wrapper(F) if AA.is_constant(g)]
    # G = CGBMain(F, S, red)
    A = AA.base_ring(AA.parent(F[1]))
    G = CGBMain(F, empty([one(A)]), red)
    return unique(vcat([G_ for (_, _, G_) ∈ G]...))
end

function CGS_faithful(F::Vector{RE}, red=true) where {RE<:AA.MPolyRingElem}
    A = AA.base_ring(AA.parent(F[1]))
    G = CGBMain(F, empty([one(A)]), red)
    return unique(G)
end


#TODO: make this preserve the monomial order of F and S
## Figured it out: use @polyvar x y monomial_order=order
## However, due to a badly written macro in DynamicPolynomials,
## `order` gets expanded to DynamicPolynomials.order
## Probably need to file a pull request to get that fixed
#TODO: make this work with AbstractAlgebra.jl polynomials
function CGBMain(F::Vector{AP}, S::Vector{AP}, params::Tuple, ordering=nothing) where {AP<:AbstractPolynomial}
    if 1 ∈ S || (!isempty(S) && 1 ∈ groebner(S))
        return empty([(F, F, F)])
    else
        vars = variables(F)
        @polyvar __t __x[1:length(vars)] # Bem. __t > __x[1] > ...
        __params = tuple((v for v ∈ __x[end-length(params)+1:end])...)
        __vars = __x[1:end-length(params)]
        F_ = [subs(f, (v => x for (v, x) ∈ zip(vars, __x))...) for f ∈ F]
        S_ = [subs(s, (v => x for (v, x) ∈ zip(vars, __x))...) for s ∈ S]

        σ₁(g) = subs(g, __t => 1, (__xi => var_i for (__xi, var_i) ∈ zip(__x, vars))...)

        # if ordering ==
        ordering = ProductOrdering(DegLex(__vars), DegLex(collect(__params)))
        # end

        G = groebner([(__t.*F_) ; (1-__t).*S_], ordering=ProductOrdering(Lex(__t), ordering))
        hs = unique([leading_coefficient(g, __params) for g ∈ G if
              divides(__t, leading_term(g)) &&
              !has_only(leading_coefficient(g, (Tuple(__x))), __params)])
        h = reduce(lcm, hs, init=one(G[1]))
        # println(G)
        if isempty(S)
            return vcat([(S, [σ₁(h)], σ₁.(G))], [CGBMain(F, [S ; [σ₁(hi)]], params) for hi ∈ hs]...)
        else
            return vcat([(groebner(S), [σ₁(h)], σ₁.(G))], [CGBMain(F, [S ; [σ₁(hi)]], params) for hi ∈ hs]...)
        end

        # The following doesn't work, since t may have the same gensym'ed id as an existing variable
        # Consider filing a PR to DynamicPolynomials to have the variable name included in the hash
        # to avoid this problem.
        #
        # @polyvar t

        # σ₁(g) = subs(g, t => 1)

        # if ordering == nothing
        #     ordering = Lex(vars...)
        # end

        # G = groebner([(t .* F) ; (1-t) .* S], ordering=ProductOrdering(Lex(t), ordering))
        # hs = [leading_coefficient(g, params) for g ∈ G if
        #       divides(t, leading_term(g)) &&
        #       !has_only(leading_coefficient(g, Tuple(vars)), params)]
        # h = reduce(lcm, hs, init=one(G[1]))
        # return vcat([(S, [h], σ₁.(G))], [CGBMain(F, [S ; [hi]], params) for hi ∈ hs]...)

        #

    end
end

export CGBMain

function CGB(F::Vector{AP}, params::Tuple; ordering=nothing) where {AP<:AbstractPolynomial}
    # Maybe this should be a product order instead of Lex, to minimize the assumptions on the
    # original order of the variables
    S = [g for g ∈ groebner(F, ordering=Lex()) if has_only(g, params)]
    G = CGBMain(F, S, params, ordering)
    return collect(Set(vcat(S, [G_ for (_, _, G_) ∈ G]...)))
end

function CGS_faithful(F, params; ordering=nothing)
    S = [g for g ∈ groebner(F, ordering=Lex()) if has_only(g, params)]
    G = CGBMain(F, S, params, ordering)
    if isempty(S)
        return G
    else
        return [[(empty(S), S, S)] ; G]
    end
end

export CGB, CGS_faithful



function stratification(G, params, grb=false)
    if !grb
        G = CGB(G, params)
    end

    A = empty(G)

    while any([!has_only(g, params) && !isconstant(leading_coefficient(g, params)) for g ∈ G])
        println(G)
        H = [leading_coefficient(g, params) for g in G
                 if !has_only(g, params) &&
                     !isconstant(leading_coefficient(g, params))
                     ]
        h = groebner(H)[1]
        G = groebner([G ; [h]], ordering=Lex())
        A = [A ; [h]]
    end

    return A

    # println(A)

    # A_ = [A[1]]
    # for f ∈ A
    #     if groebner(A_) != groebner([A_ ; [f]])
    #         A_ = [A_ ; [f]]
    #     end
    # end
    # return A_
end

export stratification
