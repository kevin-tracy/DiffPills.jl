using LinearAlgebra
using StaticArrays
import ForwardDiff as fd
using Convex
using ECOS
let


    a = randn(3)
    b = randn(3)
    c = randn(3)

    d = Variable(3)
    θ = Variable()

    prob = minimize(norm(d-c))
    prob.constraints += d == θ*a + (1-θ)*b

    solve!(prob, ECOS.Optimizer)

    d1 = vec(d.value)
    θ1 = θ.value

    F = a - b
    g = b - c

    Q = 2*F'*F
    q = 2*F'*g

    θ2 = -Q\q

    @show θ1
    @show θ2

    d2 = θ2*a + (1-θ2)*b

    @show d1
    @show d2

    n = normalize(b - a)
    d3 = a + dot(c - a,n)*n
    @show d3

    n1 = a - b
    n2 = b - c
    θ4 = -dot(n1, n2) / dot(n1, n1)

    @show θ4

end
