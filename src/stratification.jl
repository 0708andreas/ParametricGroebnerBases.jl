import Oscar
import Base.show

struct GrobnerCover
    LMs
    basis
    P_rep
end

struct ParSet
    LMs
    basis
    P_rep
end

struct P_representation
    zero
    nonzero
end


function grobcov(F)
    Sing = Oscar.Singular
    cov = Sing.LibGrobcov.grobcov(F)

    [ParSet(s[1], s[2], [P_representation(seg[1], seg[2]) for seg in s[3]]) for s ∈ cov]
end



function grobcov(F::Vector{AP}, params) where {AP<:AbstractPolynomial}
    vars = [x for x ∈ variables(F) if !(x ∈ params)]
    Sing = Oscar.Singular
    R, U = Sing.FunctionField(Sing.QQ, [p.name for p ∈ params])
    S, X = Sing.polynomial_ring(R, [x.name for x ∈ vars])


end

function Base.show(io::IO, p::P_representation)
    print(io, "V(", join(Oscar.Singular.gens(p.zero), ", "), ")  \\  (",
          join(["V(" * join(Oscar.Singular.gens(s), ", ") * ")" for s in p.nonzero], " ∪ "), ")")
end


function Base.show(io::IO, p::ParSet)
    println(io, "Parametric Set.")
    print(io, "Leading monomials: ")
    println(io, '[', join(repr.(Oscar.Singular.gens(p.LMs)), ", "), ']')
    print(io, "Basis: ")
    println(io, "[", join(repr.(Oscar.Singular.gens(p.basis)), ", "), "]")
    println(io, "P-representation:")
    for (i, s) in enumerate(p.P_rep)
        println(io, "Segment ", i)
        println(io, "\t", s)
    end
end



export grobcov

function MP_to_AS(F::Vector{AP}) where {T, AP<:AbstractPolynomial{T}}
    vars = variables(F)
    if T == Int64
        R, __x = Oscar.polynomial_ring(Oscar.QQ, [v.name for v ∈ vars])
    elseif T == Rational{Int64}
        R, __x = Oscar.polynomial_ring(Oscar.QQ, [v.name for v ∈ vars])
    end
    F_ = [f((v => x for (v, x) ∈ zip(vars, __x))...) for f ∈ F]
    return (F_, R, __x)
end

function AS_to_MP(f, vars, c_type)
    cs = c_type.(Oscar.coefficients(f))
    vs = [Oscar.exponent_vector(f, i) for i ∈ 1:length(cs)]
    return sum([c * prod(vars .^ v) for (c, v) ∈ zip(cs, vs)])
end


function radical(I::Vector{AP}) where {T, AP<:AbstractPolynomial{T}}
    if isempty(I)
        return I
    end
    vars = variables(I)
    I_, R, __x = MP_to_AS(I)
    I_ = Oscar.radical(Oscar.ideal(R, I_))
    return [AS_to_MP(f, vars, T) for f ∈ Oscar.gens(I_)]
end

function strat(F, params)
    GS = [(radical(E), radical(N), G,
           (if (G == [0]); G; else groebner([polynomial(leading_monomial(g, params)) for g ∈ G]) end))
          for (E, N, G) ∈ CGS_alt2(F, params)]
    GS = [(E, N, G, T) for (E, N, G, T) ∈ GS if !Oscar.is_subset(N, E)]

    return unique(GS)
end

export strat, MP_to_AS, radical
