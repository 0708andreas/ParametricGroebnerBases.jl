# ParametricGrobnerBases.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://0708andreas.github.io/ParametricGroebnerBases.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://0708andreas.github.io/ParametricGroebnerBases.jl/dev/)
[![Build Status](https://github.com/0708andreas/ParametricGroebnerBases.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/0708andreas/ParametricGroebnerBases.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/0708andreas/ParametricGroebnerBases.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/0708andreas/ParametricGroebnerBases.jl)

This project implements the Suzuki-Sato algorithm for computing Gröbner systems and comprehensive/parametric Gröbner bases. For an introduction to its use, see the documentation. For a mathematical introduction to parametric Gröbner bases, see my [masters project](https://github.com/0708andreas/speciale).

## What is this all about? (for non-mathematicians)
Suppose you have a set of polynomial equations $f_1(x_1, x_2, ..., x_n) = 0, ..., f_k(x_1, x_2, ..., x_n) = 0$. A Gröbner basis is a standard tool, which enables you to answer questions about these equations. One such question could be "if I assume that these equations are satisfied, does it follow that this new equation, $g(x_1, x_2, ..., x_n) = 0$ is satisfied?" Gröbner bases are a well-established tool to answer this question using the so-called _multivariate division algorithm_, and there are many fast packages that computes Gröbner bases. In Julia, I recommend Groebner.jl.

Sometimes, these equations has parameters, which are not known ahead of time. This throws a wrench in the usual Gröbner basis machinery. However, this package can give you a _Gröbner system_. A Gröbner system is a set of Gröbner bases, and information telling you which Gröbner basis works for which choices of parameters. This enables you to answer the same kinds of questions as above using the _multivariate pseudo-division algorithm_, which is also implemented in this package.

## What is this all about? (for mathematicians)
Let $I \subset K[U][X]$ be a multivariate polynomial ideal in some variables $X$ and parameters $U$. Let $\sigma : K[U] \to K$ be a ring homomorphism such that $\sigma|_K = id$. We wish to study the behaviour of the reduced Gröbner basis of $\langle \sigma(I) \rangle$ for all ring homs $\sigma$. This package computes a _Grôbner system_, which is a finite set of triples $(E, N, G)$. These triples satisfy that $\sigma(G)$ is the reduced Gröbner basis of $\langle \sigma(I) \rangle$ if $\sigma(e) = 0$ for all $e \in E$ and $\sigma(n) \neq 0$ for all $n \in N$.

This package also implements _pseudo-division_. Pseudo-division is akin to normal multivariate division, except we allow scaling the polynomial. For example, pseudo-reducing $ax + 1$ by $bx + y$ is $b(ax + 1) = a(bx + y) - ay + b$. Hence, the pseudo-remainder is $-ay + 1$. For a precise definition of pseudo-division and its properties, see my [masters project](https://github.com/0708andreas/speciale).
