using Groebner
import MultivariatePolynomials.leading_coefficient,
    MultivariatePolynomials.ordering,
    MultivariatePolynomials.leading_monomial
import Base.lcm

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

# function CGSMain(F::Vector{AP}, params) where {AP<:AbstractPolynomial}
#     G = groebner(F)
#     @assert isgroebner(G)
#     if 1 ∈ G
#         return [(one(F[1]), F)] # Skal det være F eller G til sidst?
#     else
#         hs = [leading_coefficient(g, params) for g ∈ G
#               if !has_only(g, params)]
#         h = reduce(lcm, hs)
#         # if h == 1
#         #     return [(h, G)]
#         # else
#             return vcat([(h, G)], [CGSMain([G ; [hi]], params) for hi in hs]...)
#         # end
#     end
# end

function CGSMain(F, params)
    G = groebner(F, ordering=Lex())
    hs = [leading_coefficient(g, params) for g ∈ G
              if !has_only(g, params)]
    h = reduce(lcm, hs, init=one(F[1]))
    if h == 1
        return [(h, F)]
    else
        return vcat([(h, G)], [CGSMain([G ; [hi]], params) for hi in hs]...)
    end
end


function CGS(F::Vector{AP}, params::Tuple) where {AP<:AbstractPolynomial}
    G₀ = groebner(F, ordering=Lex())
    H = CGSMain(G₀, params)
    GS = empty([(F, F, F)])
    if any([has_only(g, params) for g ∈ G₀])
        GS = [(empty(F), filter(g -> has_only(g, params), G₀), [one(F[1])])]
    end
    for (h, G) ∈ H
        if !any(isconstant.(filter(g -> has_only(g, params), G)))
            push!(GS, (groebner(filter(g -> has_only(g, params), G)), [h], (filter(g -> !has_only(g, params), G))))
        end
    end
    # @assert all(isgroebner.([G for (_, _, G) ∈ GS]))
    return GS
end

function CGS_alt(F::Vector{AP}, params::Tuple) where {AP<:AbstractPolynomial}
    function CGS_inner(F_, E)
        if any(isconstant.(filter(g -> has_only(g, params), E)))
            return empty([(E, E, E)])
        end
        G = groebner(F_, ordering=Lex())
        @assert isgroebner(G, ordering=Lex())
        hs = [leading_coefficient(g, params) for g ∈ G
                  if !has_only(g, params)]
        h = reduce(lcm, hs, init=one(F[1]))
        # if h == 1
            # return [(E, [h], G)]
        # else
            return vcat([(E, [h], G)],
                        [(CGS_inner([G ; [hi]], groebner([E ; [hi]]))) for hi ∈ hs]...)
        # end
    end
    return CGS_inner(F, empty(F))
end

function CGS_alt2(F, params)
    function CGS_inner(S)
        if 1 ∈ S || (!isempty(S) && 1 ∈ groebner(S))
            return empty([(F, F, F)])
        end

        G = groebner([F ; S], ordering=Lex())

        @assert isgroebner(G, ordering=Lex())
        hs = [leading_coefficient(g, params) for g ∈ G
                  if !has_only(g, params)]
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

    S = [g for g ∈ groebner(F, ordering=Lex()) if has_only(g, params)]
    if isempty(S)
        return CGS_inner(S)
    else
        return [[(empty(F), S, [one(F[1])])] ; CGS_inner(S)]
    end
end


export CGS_alt, CGS_alt2

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
        F_ = [subs(f, (v => x for (v, x) ∈ zip(vars, __x))...) for f ∈ F]
        S_ = [subs(s, (v => x for (v, x) ∈ zip(vars, __x))...) for s ∈ S]

        σ₁(g) = subs(g, __t => 1, (__xi => var_i for (__xi, var_i) ∈ zip(__x, vars))...)

        # if ordering == nothing
            ordering = Lex(__x)
        # end

        G = groebner([(__t.*F_) ; (1-__t).*S_], ordering=ProductOrdering(Lex(__t), ordering))
        hs = [leading_coefficient(g, __params) for g ∈ G if
              divides(__t, leading_term(g)) &&
              !has_only(leading_coefficient(g, (Tuple(__x))), __params)]
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
    return [[(empty(S), S, S)] ; G]
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
