

# check attitude stuff
for i = 1:1000
    p = @SVector randn(3)

    Q1 = DP.dcm_from_p(p)

    sp = SA[0 -p[3] p[2]; p[3] 0 -p[1]; -p[2] p[1] 0]
    Q2 = I + (8*sp^2 + 4*(1 - dot(p,p))*sp)/(1 + dot(p,p))^2

    @test norm(Q1-Q2) < 1e-10
end

# check dcm jacobian
for i = 1:1000
    p = @SVector randn(3)
    v = @SVector randn(3)

    Q = DP.dcm_from_p(p)

    J1 = DP.mrp_dcm_jacobian(p,v)
    J2 = FD.jacobian(_p -> DP.dcm_from_p(_p)*v, p)

    @test norm(J1 - J2) < 1e-12
end

# check dcm transpose jacobian
for i = 1:1000
    p = @SVector randn(3)
    v = @SVector randn(3)

    Q = DP.dcm_from_p(p)

    J1 = DP.mrp_dcm_T_jacobian(p,v)
    J2 = FD.jacobian(_p -> DP.dcm_from_p(_p)'*v, p)

    @test norm(J1 - J2) < 1e-12
end

# check dcm tilde jacobian
for i = 1:1000
    p = @SVector randn(3)
    v = @SVector randn(2)

    J1 = DP.mrp_dcm_tilde_jacobian(p,v)
    J2 = FD.jacobian(_p -> DP.dcm_tilde_from_p(_p)*v, p)

    @test norm(J1 - J2) < 1e-12
end

# check dcm tilde transpose jacobian
for i = 1:1000
    p = @SVector randn(3)
    v = @SVector randn(3)

    J1 = DP.mrp_dcm_tilde_T_jacobian(p,v)
    J2 = FD.jacobian(_p -> DP.dcm_tilde_from_p(_p)'*v, p)

    @test norm(J1 - J2) < 1e-12
end

# check dcm tilde
for i = 1:1000
    p = @SVector randn(3)

    Q1 = DP.dcm_from_p(p)[:,1:2]

    Q2 = DP.dcm_tilde_from_p(p)

    @test norm(Q1 - Q2) < 1e-10
end

# dcm product vec(Q1'*Q2) wrt p1
for i = 1:1000
    p1 = @SVector randn(3)
    p2 = @SVector randn(3)

    Q1 = DP.dcm_from_p(p1)
    Q2 = DP.dcm_from_p(p2)

    J1 = FD.jacobian(_p1 -> vec(DP.dcm_from_p(_p1)'*Q2), p1)
    J2 = DP.dcm_prod_jac_p1(p1,Q2)

    @test norm(J1-J2) < 1e-10
end

# dcm product vec(Q1'*Q2) wrt p2
for i = 1:1000
    p1 = @SVector randn(3)
    p2 = @SVector randn(3)

    Q1 = DP.dcm_from_p(p1)
    Q2 = DP.dcm_from_p(p2)

    J1 = FD.jacobian(_p2 -> vec(Q1'*DP.dcm_from_p(_p2)), p2)
    J2 = DP.dcm_prod_jac_p2(Q1,p2)

    @test norm(J1-J2) < 1e-10
end

# dcm product vec(Q̃1'*Q̃2) wrt p1
for i = 1:1000
    p1 = @SVector randn(3)
    p2 = @SVector randn(3)

    Q1 = DP.dcm_tilde_from_p(p1)
    Q2 = DP.dcm_tilde_from_p(p2)

    J1 = FD.jacobian(_p1 -> vec(DP.dcm_tilde_from_p(_p1)'*Q2), p1)
    J2 = DP.dcm_tilde_prod_jac_p1(p1,Q2)

    @test norm(J1-J2) < 1e-10
end

# dcm product vec(Q̃1'*Q̃2) wrt p2
for i = 1:1000
    p1 = @SVector randn(3)
    p2 = @SVector randn(3)

    Q1 = DP.dcm_tilde_from_p(p1)
    Q2 = DP.dcm_tilde_from_p(p2)

    J1 = FD.jacobian(_p2 -> vec(Q1'*DP.dcm_tilde_from_p(_p2)), p2)
    J2 = DP.dcm_tilde_prod_jac_p2(Q1,p2)

    @test norm(J1-J2) < 1e-10
end


# quaternion stuff
for i = 1:1000
    p = @SVector randn(3)
    q = DP.q_from_p(p)

    p2 = DP.p_from_q(q)

    # compare the two MRP's
    Q1 = DP.dcm_from_p(p)
    Q2 = DP.dcm_from_p(p2)

    @test norm(Q1'*Q2 - I) < 1e-12

    # compare MRP and quat
    Q1 = DP.dcm_from_p(p)
    Q2 = DP.dcm_from_q(q)

    @test norm(Q1'*Q2 - I) < 1e-12

    # check the dcm dtildes
    Q1 = DP.dcm_tilde_from_p(p)
    Q2 = DP.dcm_tilde_from_q(q)

    @test norm(Q1 - Q2) < 1e-12

end
