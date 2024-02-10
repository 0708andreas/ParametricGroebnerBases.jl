using ParametricGrobnerBases
using MultivariatePolynomials
using Test

@testset "ParametricGrobnerBases.jl" begin
    # Write your tests here.
    # p = x*y - u + v*x^2*y^3
    # p2 = x*y*u + x*y*v
    # p3 = u + v*x*y
    # p_ = make_parameters(p, (u, v))
    # p2_ = make_parameters(p2, (u, v))
    # p3_ = make_parameters(p3, (u, v))
    # c = ParametricGrobnerBases.Condition([polynomial(u)], [polynomial(v)])
    # c2 = ParametricGrobnerBases.Condition([polynomial(v)], [polynomial(u)])

    # @test color(terms(p)[1], c) == green
    # @test color(terms(p)[2], c) == white
    # @test color(terms(p)[3], c) == red

    # @test is_reducible(c2, p_, [p2_, p3_])[2:3] == ((-1.0), true)
    # @test is_reducible(c, p_, [p3_, p2_])[2:3] == (x*y^2, true)
    # @test is_reducible(c, p_, [p2_, p3_])[2:3] == (x^2 * y^3, true)

    # @test SPol(c2, p_, p2_) == make_parameters((-v - u^2) + (u*v*)*x^2 * y^3, (u, v))
    #

    @test length(CGS_alt([x^3 - u, y^4 - v, x + y], (u, v))) == 8
    @test CGS_alt([1 + y^2, y + u*x^2], (u,)) == [([], [u], [1+y^2, y+x^2*u]), ([u], [1], [1])]
end
