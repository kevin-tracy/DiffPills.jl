
function fd_prox_capsule(P1::DP.Capsule{T},P2::DP.Capsule{T},r1,r2,p1,p2) where {T}
    P1.r = r1
    P2.r = r2
    P1.p = p1
    P2.p = p2
    [DP.proximity(P1,P2; pdip_tol = 1e-12)]
end



for i = 1:1

    # create the capsule objects
    P1 = DP.create_capsule(:MRP)
    P2 = DP.create_capsule(:MRP)

    P2.r = SA[0,3.0,0]
    P1.R = 0.75;                P1.R = 1.5
    P1.L = 3.0;                 P1.L = 2.0

    # calculate proximity
    a,b = DP.get_ends(P1)
    c,d = DP.get_ends(P2)

    @test DP.check_parallel(a,b,c,d)

    Q, q, G, h = DP.get_cost_terms(a,b,c,d)

    z1 = DP.active_set_qp(Q,q;get_duals = false)
    z2,_,_ = DP.pdip(Q,q,G,h)

    @test sum(isnan.(z1)) == 2
    @test sum(isnan.(z2)) == 0

    ℓ = DP.proximity(P1, P2)

    @test sum(isnan.(ℓ)) == 0
end

for i = 1:1000
    # create the capsule objects
    P1 = DP.create_capsule(:MRP)
    P2 = DP.create_capsule(:MRP)

    # P2.r = SA[0,3.0,0]
    P1.R = 0.75;                P2.R = 1.5
    P1.L = 3.0;                 P2.L = 2.0

    r1 = 2*(@SVector randn(3))
    r2= 2*(@SVector randn(3))
    p1 = (@SVector randn(3))
    p2 = (@SVector randn(3))
    P1.r = r1
    P2.r = r2
    P1.p = p1
    P2.p = p2

    # calculate proximity gradients
    dℓ_dr1, dℓ_dp1, dℓ_dr2, dℓ_dp2 = DP.proximity_jacobians(P1,P2)
    # @show dℓ_dr1

    fd_dℓ_dr1 = FiniteDiff.finite_difference_jacobian(_r1 -> fd_prox_capsule(P1,P2,_r1,r2,p1,p2), r1)
    fd_dℓ_dr2 = FiniteDiff.finite_difference_jacobian(_r2 -> fd_prox_capsule(P1,P2,r1,_r2,p1,p2), r2)
    fd_dℓ_dp1 = FiniteDiff.finite_difference_jacobian(_p1 -> fd_prox_capsule(P1,P2,r1,r2,_p1,p2), p1)
    fd_dℓ_dp2 = FiniteDiff.finite_difference_jacobian(_p2 -> fd_prox_capsule(P1,P2,r1,r2,p1,_p2), p2)

    @test size(dℓ_dr1) == size(fd_dℓ_dr1)
    @test size(dℓ_dp1) == size(fd_dℓ_dp1)
    @test size(dℓ_dr2) == size(fd_dℓ_dr2)
    @test size(dℓ_dp2) == size(fd_dℓ_dp2)

    @test isapprox(dℓ_dr1, fd_dℓ_dr1, atol = 1e-3)
    @test isapprox(dℓ_dp1, fd_dℓ_dp1, atol = 1e-3)
    @test isapprox(dℓ_dr2, fd_dℓ_dr2, atol = 1e-3)
    @test isapprox(dℓ_dp2, fd_dℓ_dp2, atol = 1e-3)

end
