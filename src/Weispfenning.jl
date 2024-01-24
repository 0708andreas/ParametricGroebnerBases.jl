@enum Color green=0 red=1 white=2

struct Condition{AP1<:AbstractPolynomial, AP2<:AbstractPolynomial}
    rd :: Vector{AP1}
    gr :: Vector{AP2}
end

Cover = Vector{Condition}

function Base.show(io::IO, c::Condition)
    if __show_poly_types
        println(io, "Green: $(c.gr)")
        println(io, "Red:   $(c.rd)")
    else
        print(io, "(rd: [$(join(c.rd, ','))], gr: [$(join(c.gr, ','))])")
    end
end

function Base.show(io::IO, cs::Vector{C}) where {C<:Condition}
    print(io, "[$(join(cs, ','))]")
end

function Base.show(io::IO, t::Tuple{C, Vector{AP}}) where {C<:Condition, AP<:AbstractPolynomial}
    if __show_poly_types
        print(io, "($(t[1]), $(dump(t[2])))")
    else
        print(io, "($(t[1]), [$(join(t[2], ','))])")
    end
end



function color(p, cond::Condition)
    if iszero(p)
        return green
    elseif maxdegree(p) == 0 # p is a non-zero constant polynomial
        return red
    elseif iszero(rem(p, cond.gr)) || any([iszero(rem(p, p_)) for p_ in cond.gr])
        return green
    # elseif any([iszero(rem(p, p_)) for p_ in cond.rd])
    elseif any([(maxdegree(rem(p, p_)) == 0) && (maxdegree(div(p, p_)) == 0) for p_ in cond.rd])
        return red
    else
        return white
    end
end

function MP.leading_term(p::AP, γ::Condition) where {AP<:AbstractPolynomial}
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
        while ((t, c) = leading_term(p, γ))[2] == white
            push!(Γ_, Condition([γ.rd ; [coefficient(t)]], γ.gr))
            γ = Condition(γ.rd, [γ.gr ; [coefficient(t)]])
        end
        push!(Γ_, γ)
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
    coef = div(coefficient(t1), coefficient(t2))
    mon = monomial(MP.leading_term(div(monomial(t1), monomial(t2))))
    return term(coef, mon)
    # term(div(coefficient(t1), coefficient(t2)), div(monomial(t1), monomial(t2)))
end

function Base.lcm(a::T1, b::T2) where {AP1<:AbstractPolynomial, AP2<:AbstractPolynomial,
                                       T1<:AbstractTerm{AP1}, T2<:AbstractTerm{AP2}}
    coef = lcm(coefficient(a), coefficient(b))
    mon = lcm(monomial(a), monomial(b))
    return term(coef, mon)
end

function PART(γ, p::AP) where {AP<:AbstractPolynomial}
    p_ = zero(p)
    for t ∈ terms(p)
        if color(t, γ) != green
            p_ = p_ + t
        end
    end
    return p_
end

function SPol(γ::Condition, f::AP, g::AP) where {AP<:AbstractPolynomial}
    t_f, c_f = leading_term(f, γ)
    t_g, c_g = leading_term(g, γ)
    @assert c_f == red "c_f = $c_f"
    @assert c_g == red "c_g = $c_g"

    a, s = coefficient(t_f), monomial(t_f)
    b, t = coefficient(t_g), monomial(t_g)
    w = lcm(s, t)
    d = div(w, s)
    e = div(w, t)

    return map_coefficients(p -> b*p, d*f) - map_coefficients(p -> a*p, e*g)
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
        println("----------------")
        println(length(P))
        γ, G, f, g = pop!(P)
        t_f, c_f = leading_term(f, γ)
        t_g, c_g = leading_term(g, γ)
        if c_f == green || c_g == green
            continue
        end
        delete!(GS, (γ, G))
        println("G = $(join(G, ','))")
        println("γ = $γ")
        println("f = $(PART(γ, f))")
        println("g = $(PART(γ, g))")
        h = SPol(γ, f, g)
        println("h = $(PART(γ, h))")
        k, c = NORMALFORM(γ, h, G)
        println("k = $(PART(γ, k))")
        Δ = DET(γ, [k])
        println("Δ = $Δ")

        Δ_ = [δ for δ ∈ Δ if has_red_terms(k, δ)]
        Δ_C = [δ for δ ∈ Δ if !has_red_terms(k, δ)]
        println("Δ_ = $(Δ_)")
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
