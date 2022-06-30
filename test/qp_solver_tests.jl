
seed!(1234)


for i = 1:10000
    F = @SMatrix randn(3,2)

    Q = F'*F

    q = @SVector randn(2)

    G = SA[1 0; 0 1; -1 0; 0 -1.0]
    h = SA[1,1,0,0.0]

    x1,_,λ1 = DP.pdip(Q,q,G,h; verbose = false, tol = 1e-14)
    x2,λ2 = DP.active_set_qp(Q, q; get_duals = true)

    @test norm(x1 - x2) < 1e-8
    @test norm(λ1 - λ2) < 1e-8
end

# test optnet
function ℓ_fd_optnet(Qv,q,G,h, dℓ_dz)
    # reshape Q
    # NOTE: when I finite diff this way, the answer is correct, but when
    # I assemble it the following way, things are off by a factor of 2
    # Q = SA[Qv[1] Qv[2];
           # Qv[2] Qv[3]]
    # could the optnet paper be wrong?
    
    Q = reshape(Qv,2,2)
    Q = 0.5*(Q' + Q)
    Q = SMatrix{2,2}(Q)
    z,_,_ = DP.pdip(Q,q,G,h; tol = 1e-12)
    SA[dℓ_dz*z]
end

for i = 1:1000
    F = @SMatrix randn(2,2)
    Q = F'*F
    Qv = vec(Q)
    q = @SVector randn(2)
    G = SA[1 0; 0 1; -1 0; 0 -1.0]
    h = SA[100,100,0,0.0]
    z,s,λ = DP.pdip(Q,q,G,h; tol = 1e-10)
    dℓ_dz = SA[1,4.0]'

    # finite diff
    JQv = FiniteDiff.finite_difference_jacobian(_Qv -> ℓ_fd_optnet(_Qv,q,G,h,dℓ_dz), Qv)
    Jq = FiniteDiff.finite_difference_jacobian(_q -> ℓ_fd_optnet(Qv,_q,G,h,dℓ_dz), q)

    # optnet
    dℓ_dq, dℓ_dQv = DP.fast_optnet(Q,G,z,s,λ,vec(dℓ_dz))

    @test size(JQv) == size(dℓ_dQv)
    @test size(Jq) == size(dℓ_dq)

    if norm(Jq)>1e-2
        @test isapprox(vec(Jq),vec(dℓ_dq), rtol = 1e-3)
    end

    if norm(JQv)>1e-2
        @test isapprox(vec(JQv),vec(dℓ_dQv), rtol = 1e-3)
    end


end
