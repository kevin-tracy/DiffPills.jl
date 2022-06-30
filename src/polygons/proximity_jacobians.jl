
@inline function Qv_jacobians(p1::SVector{3,T},
                              p2::SVector{3,T},
                              Q̃1::SMatrix{3,2,T,6},
                              Q̃2::SMatrix{3,2,T,6}) where {T}

    # jacobians of vec(Q) wrt p1 and p2

    # jacobian of Q1'*Q2 wrt p1
    Jabcd = dcm_tilde_prod_jac_p1(p1,Q̃2)
    zr = SA[0 0 0.0]
    dQv_dp1 = [
            zr;
            zr;
            Jabcd[1,:]';
            Jabcd[3,:]';
            zr;
            zr;
            Jabcd[2,:]';
            Jabcd[4,:]';
            Jabcd[1,:]';
            Jabcd[2,:]';
            zr;
            zr;
            Jabcd[3,:]';
            Jabcd[4,:]';
            zr;
            zr
        ]

        Jabcd = dcm_tilde_prod_jac_p2(Q̃1,p2)

        dQv_dp2 = [
        zr;
        zr;
        Jabcd[1,:]';
        Jabcd[3,:]';
        zr;
        zr;
        Jabcd[2,:]';
        Jabcd[4,:]';
        Jabcd[1,:]';
        Jabcd[2,:]';
        zr;
        zr;
        Jabcd[3,:]';
        Jabcd[4,:]';
        zr;
        zr
        ]

        return dQv_dp1, dQv_dp2
end

@inline function q_jacobians(p1::SVector{3,T},
                             p2::SVector{3,T},
                             r1::SVector{3,T},
                             r2::SVector{3,T},
                             Q̃1::SMatrix{3,2,T,6},
                             Q̃2::SMatrix{3,2,T,6}) where {T}
    dq_dp1 = vcat(mrp_dcm_tilde_T_jacobian(p1,r1-r2),
                  SA[0 0 0.0; 0 0 0.0])
    dq_dp2 = vcat(SA[0 0 0.0; 0 0 0.0],
                  -mrp_dcm_tilde_T_jacobian(p2,r1-r2))

    dq_dr1 = vcat(Q̃1', -Q̃2')
    dq_dr2 = -dq_dr1

    return dq_dp1, dq_dp2, dq_dr1, dq_dr2
end

function poly_ℓ_derivs(r1,r2,R1,R2,p1,p2,Q̃1, Q̃2, z)

    z1 = SA[z[1],z[2]]
    z2 = SA[z[3],z[4]]
    e1 = r1 + Q̃1*z1
    e2 = r2 + Q̃2*z2
    e = e1-e2

    dℓ_dz = -[-2*Q̃1'*e;2*Q̃2'*e]'

    dℓ_dr1 = 2*e'
    dℓ_dr2 = -2*e'

    # dℓ_drot1 =  2*e'
    # dℓ_drot2 = -2*e'

    drot1_dp1 = mrp_dcm_tilde_jacobian(p1,z1)
    drot2_dp2 = mrp_dcm_tilde_jacobian(p2,z2)

    dℓ_dp1 = dℓ_dr1*drot1_dp1 # replaced dℓ_drot1 with dℓ_dr1
    dℓ_dp2 = dℓ_dr2*drot2_dp2 # replaced dℓ_drot2 with dℓ_dr2

    return dℓ_dz, dℓ_dr1, dℓ_dr2, dℓ_dp1, dℓ_dp2
end


@inline function proximity_jacobians(P1::Polygon{T,nh_1,nh2_1},
                                     P2::Polygon{T,nh_2,nh2_2}; pdip_tol = 1e-9) where {T,nh_1,nh_2,nh2_1,nh2_2}

    Q̃1 = dcm_tilde_from_p(P1.p)
    Q̃2 = dcm_tilde_from_p(P2.p)
    #
    Q, q, G, h = get_qp_terms(P1,P2,Q̃1,Q̃2)
    #
    z,s,λ = pdip(Q,q,G,h; tol = pdip_tol)

    dQv_dp1, dQv_dp2 = Qv_jacobians(P1.p, P2.p, Q̃1, Q̃2)

    dq_dp1, dq_dp2, dq_dr1, dq_dr2 = q_jacobians(P1.p,P2.p,P1.r,P2.r, Q̃1, Q̃2)

    dℓ_dz, dℓ_dr1, dℓ_dr2, dℓ_dp1, dℓ_dp2 = poly_ℓ_derivs(P1.r,P2.r,P1.R,P2.R,P1.p,P2.p,Q̃1, Q̃2, z)

    # full matrix method (works for active set and pdip)
    # dℓ_dq, dℓ_dQv = optnet(z,λ,vec(dℓ_dz),Q,G,h)

    # fast version (only works for pdip)
    dℓ_dq, dℓ_dQv = fast_optnet(Q,G,z,s,λ,vec(dℓ_dz))

    dℓ_dr1 = dℓ_dr1 + dℓ_dq*dq_dr1
    dℓ_dr2 = dℓ_dr2 + dℓ_dq*dq_dr2

    dℓ_dp1 = dℓ_dp1 + dℓ_dQv*dQv_dp1 + dℓ_dq*dq_dp1
    dℓ_dp2 = dℓ_dp2 + dℓ_dQv*dQv_dp2 + dℓ_dq*dq_dp2

    return dℓ_dr1, dℓ_dp1, dℓ_dr2, dℓ_dp2
end
