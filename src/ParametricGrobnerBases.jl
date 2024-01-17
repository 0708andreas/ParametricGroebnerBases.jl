module ParametricGrobnerBases

using MultivariatePolynomials
MP = MultivariatePolynomials
using DynamicPolynomials
using GLMakie
import Contour
import Base.div

@polyvar x y u v

@enum Color green=0 red=1 white=2

struct Condition{AP1<:AbstractPolynomial, AP2<:AbstractPolynomial}
    gr :: Vector{AP1}
    rd :: Vector{AP2}
end

Cover = Vector{Condition}

function color(p, cond::Condition)
    if iszero(p)
        return green
    elseif maxdegree(p) == 0 # p is a non-zero constant polynomial
        return red
    elseif any([iszero(rem(p, p_)) for p_ in cond.gr])
        return green
    elseif any([iszero(rem(p, p_)) for p_ in cond.rd])
        return red
    else
        return white
    end
end

function leading_term(p::AP, γ::Condition) where {AP<:AbstractPolynomial}
    for t in reverse(terms(p))
        if color(coefficient(t), γ) == red
            return (t, red)
        elseif color(coefficient(t), γ) == white
            return (t, white)
        end
    end
    return (MP.leading_term(p), green)
end

function has_red_terms(p::AP, γ::Condition) where {AP<:AbstractPolynomial}
    for t in terms(p)
        if color(coefficient(t), γ) == red
            return true
        end
    end
    return false
end

function extend_cover(Γ, p::AP) where {AP<:AbstractPolynomial}
    Γ_ = empty(Γ)
    for γ ∈ Γ
        t, c = leading_term(p, γ)
        if c == white
            push!(Γ_, Condition([γ.gr ; [coefficient(t)]], γ.rd))
            push!(Γ_, Condition(γ.gr, [γ.rd ; [coefficient(t)]]))
        else
            push!(Γ_, γ)
        end
    end
    return Γ_
end

function DET(γ::Condition, F::Vector{AP}) where {AP<:AbstractPolynomial}
    Γ = [γ]
    for p in F
        Γ = extend_cover(Γ, p)
    end
    return Γ
end

# function is_reducible(γ, f, P)
#     t_f, c_f = leading_term(f, γ)
#     if c_f == green
#         return (f, false)
#     end
#     for p in P
#         t_p, c_p = leading_term(p, γ)
#         if iszero(rem(t_f, t_p))
#             return (p, true)
#         end
#     end
#     return (f, false)
# end

function divides(t1::Term{AP1}, t2::Term{AP2}) where {AP1<:AbstractPolynomial, AP2<:AbstractPolynomial}
    iszero(rem(coefficient(t2), coefficient(t1))) && iszero(rem(monomial(t2), monomial(t1)))
end

function div(t1::Term{AP1}, t2::Term{AP2}) where {AP1<:AbstractPolynomial, AP2<:AbstractPolynomial}
    coef = MP.leading_term(div(coefficient(t1), coefficient(t2)))
    mon = monomial(MP.leading_term(div(monomial(t1), monomial(t2))))
    return term(coef, mon)
    # term(div(coefficient(t1), coefficient(t2)), div(monomial(t1), monomial(t2)))
end

function SPol(γ::Condition, f::AP, g::AP) where {AP<:AbstractPolynomial}
    t_f, c_f = leading_term(f, γ)
    t_g, c_g = leading_term(g, γ)
    @assert c_f == red "c_f = $c_f"
    @assert c_g == red "c_g = $c_g"

    w = lcm(t_f, t_g)
    s = div(w, t_f)
    r = div(w, t_g)

    return s*f - r*g
end

function is_reducible(γ, f, P)
    for t_f in reverse(terms(f))
        if color(coefficient(t_f), γ) == green
            continue
        end
        for p in P
            t_p, c_p = leading_term(p, γ)
            if divides(t_p, t_f)
                return (p, div(t_f, t_p), true)
            end
        end
    end
    return (f, f, false)
end

function NORMALFORM(γ, f, P) # Assumes that γ determines P
    g = f
    c = zero(f) + 1
    while ((p, q, b) = is_reducible(γ, g, P))[3]
        t_p, c_p = leading_term(p, γ)
        g = g - q*p
        c = c*coefficient(t_p)
    end
    return (g, c)
end

function GROBNERSYSTEM(F, B)
    Γ = vcat([DET(γ, F) for γ ∈ B]...)
    GS = Set([(γ, F) for γ ∈ Γ])
    P = Set([(γ, F, f, g) for γ ∈ Γ for f ∈ F for g ∈ F if f != g])
    while !isempty(P)
        println(length(P))
        γ, G, f, g = pop!(P)
        delete!(GS, (γ, G))
        t_f, c_f = leading_term(f, γ)
        t_g, c_g = leading_term(g, γ)
        if c_f == green || c_g == green
            continue
        end
        
        h = SPol(γ, f, g)
        k, c = NORMALFORM(γ, h, G)
        Δ = DET(γ, [k])
        Δ_ = [δ for δ ∈ Δ if has_red_terms(k, δ)]
        Δ_C = [δ for δ ∈ Δ if !has_red_terms(k, δ)]
        if isempty(Δ_)
            push!(GS, (γ, G))
        else
            union!(GS, Set([(δ, [G ; [k]]) for δ ∈ Δ_]))
            union!(GS, Set([(δ, G) for δ ∈ Δ_C]))

            P_ = filter(p -> p[1] == γ && p[2] == G, P)
            filter!(p -> p[1] != γ || p[2] != G, P)
            union!(P, Set([(δ, [G ; [k]], f_, k) for f_ ∈ G for δ ∈ Δ_]))
            union!(P, Set([(δ, G, f_, g_) for (_, _, f_, g_) ∈ P_ for δ ∈ Δ]))
        end
    end
    return GS
end

function make_parameters(p, vars)
    p_ = polynomial(zero(p), polynomial_type(leading_coefficient(p)*prod(vars)))
    for t in terms(p)
        mon = monomial(subs(t, [var => 1 for var in vars]...))

        remaining_vars = filter(x -> !(x in vars), variables(monomial(t)))
        coef = subs(t, [var => 1 for var in remaining_vars]...)

        p_ = p_ + term(coef, mon)
    end
    return p_
end

unparameterize(t) = coefficient(t) * monomial(t)


function plot_poly(p, params)
    f = Figure()
    ax = Axis(f[1, 1])
    sliders = []
    p_ = subs(p, (par => 5 for par in params)...)
    # for (i, par) in enumerate(params)
    #     sl = Slider(f[1, i+1], range=-10:0.01:10, startvalue=5)
    #     on(sl.value) do _
    #         p_[] = subs(p, (par => slider.value[] for (par, slider) in zip(params, sliders))...)
    #     end

    #     push!(sliders, Slider(f[1, i+1], range=-10:0.01:10, startvalue=5))
    # end

    r = -10:0.1:10
    z = [p_(x => xi, y => yi) for xi in r, yi in r]
    lns = Contour.lines(Contour.contour(collect(r), collect(r), z, 0.0))
    empty!(ax)
    for ln in lns
        xs, ys = Contour.coordinates(ln)
        lines!(ax, xs, ys)
    end
    return f
end


export x, y, u, v, make_parameters, Color, green, red, white, Condition, color,
    leading_term, extend_cover, DET, is_reducible, NORMALFORM, SPol, GROBNERSYSTEM

end
