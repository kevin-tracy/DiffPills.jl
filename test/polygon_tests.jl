for i = 1:1000
    p1 = @SVector randn(3)
    p2 = @SVector randn(3)

    Q1 = DP.dcm_tilde_from_p(p1)
    Q2 = DP.dcm_tilde_from_p(p2)

    prod_12 = Q1'*Q2
    a,b,c,d = prod_12

    prod_21 = Q2'*Q1
    # [a c
     # b d]  is prod_12

    # [a b
     # c d]  is prod_21

    Q = [I(2) prod_12; prod_21 I(2)]
    v1 = vec(Q)
    v2 = [
            1;
            0;
            a;
            c;
            0;
            1;
            b;
            d;
            a;
            b;
            1;
            0;
            c;
            d;
            0;
            1
        ]

    @test norm(v1-v2) < 1e-12

end

for i = 1:1000
    p1 = @SVector randn(3)
    p2 = @SVector randn(3)

    Q1 = DP.dcm_tilde_from_p(p1)
    Q2 = DP.dcm_tilde_from_p(p2)

    prod_12 = Q1'*Q2
    a,b,c,d = prod_12

    # [a c
     # b d]  is prod_12

    J1 = FD.jacobian(_p1 -> vec([I(2)   DP.dcm_tilde_from_p(_p1)'*DP.dcm_tilde_from_p(p2);
         DP.dcm_tilde_from_p(p2)'*DP.dcm_tilde_from_p(_p1) I(2)]), p1)


    Jabcd = DP.dcm_tilde_prod_jac_p1(p1,Q2)
    zr = SA[0 0 0.0]
    J2 = [
            zr;
            zr;
            Jabcd[1,:]';# 1
            Jabcd[3,:]';# 3
            zr;
            zr;
            Jabcd[2,:]';# 2
            Jabcd[4,:]';# 4
            Jabcd[1,:]';# 1
            Jabcd[2,:]';# 2
            zr;
            zr;
            Jabcd[3,:]'; # 3
            Jabcd[4,:]'; # 4
            zr;
            zr
        ]
    @test norm(J1-J2) < 1e-12

end


for i = 1:1000
    p1 = @SVector randn(3)
    p2 = @SVector randn(3)

    Q1 = DP.dcm_tilde_from_p(p1)
    Q2 = DP.dcm_tilde_from_p(p2)

    prod_12 = Q1'*Q2
    a,b,c,d = prod_12

    dQv_dp1, dQv_dp2 = DP.Qv_jacobians(p1,p2,Q1,Q2)

    J1 = FD.jacobian(_p1 -> vec([I(2)   DP.dcm_tilde_from_p(_p1)'*DP.dcm_tilde_from_p(p2);
         DP.dcm_tilde_from_p(p2)'*DP.dcm_tilde_from_p(_p1) I(2)]), p1)


    @test norm(J1-dQv_dp1) < 1e-12

    J2 = FD.jacobian(_p2 -> vec([I(2)   DP.dcm_tilde_from_p(p1)'*DP.dcm_tilde_from_p(_p2);
         DP.dcm_tilde_from_p(_p2)'*DP.dcm_tilde_from_p(p1) I(2)]), p2)

    @test norm(J2-dQv_dp2) < 1e-12

end


for i = 1:1000
    p1 = @SVector randn(3)
    p2 = @SVector randn(3)

    r1 = @SVector randn(3)
    r2 = @SVector randn(3)

    Q̃1 = DP.dcm_tilde_from_p(p1)
    Q̃2 = DP.dcm_tilde_from_p(p2)

    J1 = FD.jacobian(_p1 -> [DP.dcm_tilde_from_p(_p1)'*(r1 - r2); -Q̃2'*(r1-r2)],p1)
    J2 = FD.jacobian(_p2 -> [Q̃1'*(r1 - r2); -DP.dcm_tilde_from_p(_p2)'*(r1-r2)],p2)
    J3 = FD.jacobian(_r1 -> [Q̃1'*(_r1 - r2); -Q̃2'*(_r1-r2)],r1)
    J4 = FD.jacobian(_r2 -> [Q̃1'*(r1 - _r2); -Q̃2'*(r1-_r2)],r2)

    dq_dp1, dq_dp2, dq_dr1, dq_dr2 = DP.q_jacobians(p1,p2,r1,r2, Q̃1, Q̃2)

    @test norm(J1-dq_dp1) < 1e-10
    @test norm(J2-dq_dp2) < 1e-10
    @test norm(J3-dq_dr1) < 1e-10
    @test norm(J4-dq_dr2) < 1e-10
end


# this is just for testing
function ℓ_fd_test(r1,r2,R1,R2,p1,p2,z)
    # two closest points on the two planes
    e1 = r1 + DP.dcm_tilde_from_p(p1)*SA[z[1],z[2]]
    e2 = r2 + DP.dcm_tilde_from_p(p2)*SA[z[3],z[4]]

    # vector between them
    e = e1 - e2

    # return the squared distance
    [dot(e,e) - (R1 + R2)^2]
end

for i = 1:1000

    r1 = @SVector randn(3)
    r2 = @SVector randn(3)
    p1 = @SVector randn(3)
    p2 = @SVector randn(3)
    z  = @SVector randn(4)

    R1 = abs(randn())
    R2 = abs(randn())

    Q̃1 = DP.dcm_tilde_from_p(p1)
    Q̃2 = DP.dcm_tilde_from_p(p2)

    dℓ_dz, dℓ_dr1, dℓ_dr2, dℓ_dp1, dℓ_dp2 = DP.poly_ℓ_derivs(r1,r2,R1,R2,p1,p2,Q̃1, Q̃2, z)

    @test norm(dℓ_dz - FD.jacobian(_z-> ℓ_fd_test(r1,r2,R1,R2,p1,p2,_z),z)) < 1e-12

    @test norm(dℓ_dr1 - FD.jacobian(_r1-> ℓ_fd_test(_r1,r2,R1,R2,p1,p2,z),r1)) < 1e-12
    @test norm(dℓ_dr2 - FD.jacobian(_r2-> ℓ_fd_test(r1,_r2,R1,R2,p1,p2,z),r2)) < 1e-12

    @test norm((dℓ_dp1) - FD.jacobian(_p1-> ℓ_fd_test(r1,r2,R1,R2,_p1,p2,z),p1)) < 1e-12
    @test  norm((dℓ_dp2) - FD.jacobian(_p2-> ℓ_fd_test(r1,r2,R1,R2,p1,_p2,z),p2)) < 1e-12

end
