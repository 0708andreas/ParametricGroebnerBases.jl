module ParametricGrobnerBases

using MultivariatePolynomials
MP = MultivariatePolynomials
using DynamicPolynomials
using DataStructures
# using GLMakie
import Contour
import Base.div
import Base.show

__show_poly_types = false

# @polyvar x y z u v # monomial_order=LexOrder

# include("Weispfenning.jl")

include("reduction.jl")

include("SSAlgo.jl")

export CGSMain, CGS, CGBMain, CGB, leading_monomial

# include("stratification.jl")

rand_poly(vars, density=0.5) = rand(1:10) +
    polynomial(rand(1:10, length(vars)) .* (density .> rand(length(vars))), vars) +
    polynomial(rand(1:10, length(vars), length(vars)) .* (density .> rand(length(vars), length(vars))), vars)


function final_grb(I, params)
    display(I); println(""); display(groebner(I, ordering=Lex()))
    println("")

    display(make_parameters.(groebner([f - p for (f, p) in zip(I, params)], ordering=Lex()), Ref(params)))
    I_ = subs(groebner([f - p for (f, p) in zip(I, params)], ordering=Lex()), (p => 0 for p in params)...)
    display(I_)
    display(isgroebner(I_));
    println("\n\n")
    return isgroebner(I_);
end

function sturmfels(vars, params)
    done = false
    while !done
        f = rand_poly(vars)

        done = !final_grb([f*rand_poly(vars) + 1, f], params)
    end
end



export rand_poly, final_grb, sturmfels


function make_parameters(p, vars)
    p_ = polynomial(zero(p), polynomial_type(leading_coefficient(p)*prod(vars)))
    for t in terms(p)
        mon = monomial(subs(t, [var => 1 for var in vars]...))

        remaining_vars = filter(x -> !(x in vars), variables(monomial(t)))
        if isempty(remaining_vars)
            coef = t
        else
            coef = subs(t, [var => 1 for var in remaining_vars]...)
        end

        p_ = p_ + term(coef, mon)
    end
    return p_
end

unparameterize(t::AbstractTerm) = coefficient(t) * monomial(t)
unparameterize(p::AbstractPolynomial) = sum(unparameterize.(terms(p)))


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


export x, y, z, u, v, make_parameters, Color, green, red, white, Condition, color,
    leading_term, extend_cover, DET, is_reducible, NORMALFORM, SPol, GROBNERSYSTEM

end
