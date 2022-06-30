

# check proximity function
function get_ends(x,L)
    r = SVector{3}(x[1:3])
    p = SVector{3}(x[7:9])

    n_Q_b = DP.dcm_from_p(p)
    a = r + n_Q_b*[-L/2,0,0]
    b = r + n_Q_b*[L/2,0,0]
    return a, b
end

function ineq_con_x(params,x)
    L1,L2,R1,R2 = params.L1, params.L2, params.R1, params.R2
    x1 = x[1:12]
    x2 = x[13:24]
    a,b = get_ends(x1,L1)
    c,d = get_ends(x2,L2)

    F = [(a-b) (d-c)]

    Q = F'*F
    q = F'*(b-d)

    G = [1 0; 0 1; -1 0; 0 -1]
    h = [1;1;0;0]

#     z,_,_ = mypdip(Q,q,G,h)
    z = DP.active_set_qp(Q,q;get_duals = false)

    return [DP.ℓ_fx(a,b,c,d,z,R1,R2)]
end


for i = 1:1000

    params = (
        L1 = 2 + rand(),
        R1 = 0.75 + .5*rand(),
        L2 = 1.8 + 0.8*rand(),
        R2 = 1 + 0.5*rand(),
    )

    x = @SVector randn(24) # two rigid bodies fully controlled
    x1 = x[SVector{12}(1:12)]
    r1 = SVector{3}(x1[1:3])
    p1 = SVector{3}(x1[7:9])
    x2 = x[SVector{12}(13:24)]
    r2 = SVector{3}(x2[1:3])
    p2 = SVector{3}(x2[7:9])
    q1 = DP.q_from_p(p1)
    q2 = DP.q_from_p(p2)

    pill_1 = DP.Capsule(r1,p1,q1,params.R1,params.L1,:MRP)
    pill_2 = DP.Capsule(r2,p2,q2,params.R2,params.L2,:MRP)


    c = DP.proximity(pill_1, pill_2) # scalar
    c2 = ineq_con_x(params,x)             # 1d vector (for FiniteDiff jac)

    @test norm(c-c2[1]) < 1e-12

    dℓ_dr1, dℓ_dp1, dℓ_dr2, dℓ_dp2 = DP.proximity_jacobians(pill_1,pill_2)

    dc_dx = FiniteDiff.finite_difference_jacobian(_x -> ineq_con_x(params,_x),x)
    dc_dr1 = dc_dx[:,1:3]
    dc_dp1 = dc_dx[:,7:9]
    dc_dr2 = dc_dx[:,(1:3) .+ 12]
    dc_dp2 = dc_dx[:,(7:9) .+ 12]

    @test norm(dc_dr1 - dℓ_dr1) < 1e-4
    @test norm(dc_dp1 - dℓ_dp1) < 1e-4
    @test norm(dc_dr2 - dℓ_dr2) < 1e-4
    @test norm(dc_dp2 - dℓ_dp2) < 1e-4


    # test MRP vs Quat proximity calls
    c1 = DP.proximity(pill_1, pill_2)
    pill_1.attitude = :quat
    pill_2.attitude = :quat
    c2 = DP.proximity(pill_1, pill_2)

    @test abs(c1-c2)<1e-12

end
