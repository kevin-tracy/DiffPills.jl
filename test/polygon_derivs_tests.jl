
function create_n_sided(N,d)
    ns = [ [cos(θ);sin(θ)] for θ = 0:(2*π/N):(2*π*(N-1)/N)]
    #
    A = vcat(transpose.((ns))...)
    b = d*ones(N)
    return SMatrix{N,2}(A), SVector{N}(b)
end

function fd_prox(P1,P2,r1,r2,p1,p2)
    P1.r = r1
    P2.r = r2
    P1.p = p1
    P2.p = p2
    [DP.proximity(P1,P2;pdip_tol = 1e-12)]
end



for i = 1:1000

    A1, b1 = create_n_sided(6,1.2)
    A2, b2 = create_n_sided(6,1.4)
    R1 = abs(2*randn())
    R2 = abs(2*randn())
    P1 = DP.create_polygon(A1,b1,R1,:MRP)
    P2 = DP.create_polygon(A2,b2,R2,:MRP)

    r1 = 2*(@SVector randn(3))
    r2 = 2*(@SVector randn(3))
    p1 = @SVector randn(3)
    p2 = @SVector randn(3)
    P1.r = r1
    P2.r = r2
    P1.p = p1
    P2.p = p2

    dℓ_dr1, dℓ_dp1, dℓ_dr2, dℓ_dp2 = DP.proximity_jacobians(P1,P2)

    fd_dℓ_dr1 = FiniteDiff.finite_difference_jacobian(_r1 -> fd_prox(P1,P2,_r1,r2,p1,p2), r1)
    fd_dℓ_dr2 = FiniteDiff.finite_difference_jacobian(_r2 -> fd_prox(P1,P2,r1,_r2,p1,p2), r2)
    fd_dℓ_dp1 = FiniteDiff.finite_difference_jacobian(_p1 -> fd_prox(P1,P2,r1,r2,_p1,p2), p1)
    fd_dℓ_dp2 = FiniteDiff.finite_difference_jacobian(_p2 -> fd_prox(P1,P2,r1,r2,p1,_p2), p2)

    @test size(dℓ_dr1) == size(fd_dℓ_dr1)
    @test size(dℓ_dp1) == size(fd_dℓ_dp1)
    @test size(dℓ_dr2) == size(fd_dℓ_dr2)
    @test size(dℓ_dp2) == size(fd_dℓ_dp2)

    @test isapprox(dℓ_dr1, fd_dℓ_dr1, atol = 1e-3)
    @test isapprox(dℓ_dp1, fd_dℓ_dp1, atol = 1e-3)
    @test isapprox(dℓ_dr2, fd_dℓ_dr2, atol = 1e-3)
    @test isapprox(dℓ_dp2, fd_dℓ_dp2, atol = 1e-3)


    # check quaternion vs mrp version
    ℓ1 = DP.proximity(P1,P2;pdip_tol = 1e-12)
    P1.q = DP.q_from_p(P1.p)
    P2.q = DP.q_from_p(P2.p)
    P1.attitude = :quat
    P2.attitude = :quat
    ℓ2 = DP.proximity(P1,P2;pdip_tol = 1e-12)

    @test abs(ℓ1 - ℓ2) < 1e-12
end

# # check same attitude stuff
# for i = 1:3
#
#     A1, b1 = create_n_sided(6,1.2)
#     A2, b2 = create_n_sided(6,1.4)
#     R1 = abs(2*randn())
#     R2 = abs(2*randn())
#     P1 = DP.create_polygon(A1,b1,R1,:MRP)
#     P2 = DP.create_polygon(A2,b2,R2,:MRP)
#
#     r1 = 2*(@SVector randn(3))
#     r2 = 2*(@SVector randn(3))
#     p1 = @SVector randn(3)
#     p2 = 1*p1
#     P1.r = r1
#     P2.r = r2
#     P1.p = p1
#     P2.p = p2
#     #
#     # dℓ_dr1, dℓ_dr2, dℓ_dp1, dℓ_dp2 = DP.proximity_jacobians(P1,P2)
#     #
#     # fd_dℓ_dr1 = FiniteDiff.finite_difference_jacobian(_r1 -> fd_prox(P1,P2,_r1,r2,p1,p2), r1)
#     # fd_dℓ_dr2 = FiniteDiff.finite_difference_jacobian(_r2 -> fd_prox(P1,P2,r1,_r2,p1,p2), r2)
#     # fd_dℓ_dp1 = FiniteDiff.finite_difference_jacobian(_p1 -> fd_prox(P1,P2,r1,r2,_p1,p2), p1)
#     # fd_dℓ_dp2 = FiniteDiff.finite_difference_jacobian(_p2 -> fd_prox(P1,P2,r1,r2,p1,_p2), p2)
#     #
#     # @test size(dℓ_dr1) == size(fd_dℓ_dr1)
#     # @test size(dℓ_dp1) == size(fd_dℓ_dp1)
#     # @test size(dℓ_dr2) == size(fd_dℓ_dr2)
#     # @test size(dℓ_dp2) == size(fd_dℓ_dp2)
#     #
#     # @test isapprox(dℓ_dr1, fd_dℓ_dr1, atol = 1e-3)
#     # @test isapprox(dℓ_dp1, fd_dℓ_dp1, atol = 1e-3)
#     # @test isapprox(dℓ_dr2, fd_dℓ_dr2, atol = 1e-3)
#     # @test isapprox(dℓ_dp2, fd_dℓ_dp2, atol = 1e-3)
#
#
#     # check quaternion vs mrp version
#     P1.r = r1
#     P2.r = r2
#     P1.p = p1
#     P2.p = p2
#     @show r1
#     @show r2
#     @show p1
#     @show p2
#     ℓ1 = DP.proximity(P1,P2;pdip_tol = 1e-12)
#     # P1.q = DP.q_from_p(P1.p)
#     # P2.q = DP.q_from_p(P2.p)
#     # P1.attitude = :quat
#     # P2.attitude = :quat
#     # ℓ2 = DP.proximity(P1,P2;pdip_tol = 1e-12)
#
#     # @test abs(ℓ1 - ℓ2) < 1e-12
# end
