# NOTE: moved these to utils/mrp.jl
# @inline function skew(p)
#     SA[0 -p[3] p[2]; p[3] 0 -p[1]; -p[2] p[1] 0]
# end
# @inline function mrp_dcm_jacobian(p,v)
#     # jacobian of (dcm_from_p(p)*v) wrt p
#     p1 = 8*((p'*v)*p - (p'p)*v)
#     p2 = 4*(1 - p'p)*cross(p,v)
#     p3 = (1 + p'p)^2
#     dp1 = 8*(p*v' + (p'v)*I - 2*v*p')
#     dp2 = -4*(cross(v,p)*(-2*p)' + (1 - p'p)*skew(v))
#     dp3 = 4*(1 + p'p)*p'
#     ((dp1 + dp2)*p3 - (p1 + p2)*dp3)/(p3^2)
# end
@inline function dℓ_dpoints(z::SVector{2,T},
                            a::SVector{3,T},
                            b::SVector{3,T},
                            c::SVector{3,T},
                            d::SVector{3,T}) where {T}
    z1,z2 = z
    e1 = z1*a + (1-z1)*b
    e2 = z2*c + (1-z2)*d

    dℓ_da = 2*z1*(e1-e2)'
    dℓ_db = 2*(1-z1)*(e1 - e2)'
    dℓ_dc = -2*z2*(e1-e2)'
    dℓ_dd = -2*(1-z2)*(e1 - e2)'

    t0 = e1 - e2
    ∇ℓ_z = -SA[-(2*dot(t0,a) - 2*dot(t0,b)), 2*dot(t0,c) - 2*dot(t0,d)]

    return dℓ_da, dℓ_db, dℓ_dc, dℓ_dd, ∇ℓ_z
end

@inline function dpoints_d_r_p(capsule_1::Capsule{T},
                               capsule_2::Capsule{T}) where {T}

    p1 = capsule_1.p
    L1 = capsule_1.L
    p2 = capsule_2.p
    L2 = capsule_2.L

    da_dp1 = mrp_dcm_jacobian(p1,SA[-L1/2;0;0])
    db_dp1 = mrp_dcm_jacobian(p1,SA[L1/2;0;0])
    dc_dp2 = mrp_dcm_jacobian(p2,SA[-L2/2;0;0])
    dd_dp2 = mrp_dcm_jacobian(p2,SA[L2/2;0;0])

    return da_dp1, db_dp1,  dc_dp2, dd_dp2
end

@inline function dqpcost_dx(a::SVector{3,T},
                            b::SVector{3,T},
                            c::SVector{3,T},
                            d::SVector{3,T},
                            da_dp1::SMatrix{3,3,T,9},
                            db_dp1::SMatrix{3,3,T,9},
                            dc_dp2::SMatrix{3,3,T,9},
                            dd_dp2::SMatrix{3,3,T,9}) where {T}

    zero_row = SA[0,0,0]'

    dQv_da = [2*(a-b)';(d-c)';(d-c)';zero_row ]
    dQv_db = -dQv_da
    dQv_dc = [zero_row ; (b-a)';(b-a)';2*(c-d)']
    dQv_dd = -dQv_dc


    dQv_dr1 = dQv_da + dQv_db
    dQv_dp1 = dQv_da*da_dp1 + dQv_db*db_dp1
    dQv_dr2 = dQv_dc + dQv_dd
    dQv_dp2 = dQv_dc*dc_dp2 + dQv_dd*dd_dp2

    dq_da = [(b-d)';zero_row ]
    dq_db = [(a - 2*b + d)';(d-c)']
    dq_dc = [zero_row ;(-b + d)']
    dq_dd = [(b-a)';(b-2*d+c)']


    dq_dr1 = dq_da + dq_db
    dq_dp1 = dq_da*da_dp1 + dq_db*db_dp1
    dq_dr2 = dq_dc + dq_dd
    dq_dp2 = dq_dc*dc_dp2 + dq_dd*dd_dp2

    return dQv_dr1, dQv_dp1, dQv_dr2, dQv_dp2, dq_dr1, dq_dp1, dq_dr2, dq_dp2
end

@inline function optnet(z::SVector{2,T},
                        λ::SVector{4,T},
                        ∇ℓ_z::SVector{2,T},
                        Q::SMatrix{2,2,T,4},
                        G::SMatrix{4,2,T,8},
                        h::SVector{4,T}) where {T}

    kkt_top = hcat(Q,G'*Diagonal(λ))
    # kkt_top = hcat(Q, SA[λ[1] - λ[3], λ[2] - λ[4]])
    kkt_bottom = hcat(G, Diagonal(G*z - h) )

    kkt = vcat(kkt_top,kkt_bottom)

    zero_col = @SVector zeros(length(h))
    rhs = vcat(∇ℓ_z,zero_col)

    sol = -kkt\rhs # TODO: is there not a faster way to solve this?
    # sol = -[Q G'*diagm(λ);G diagm(G*z -h)]\[∇ℓ_z;zeros(length(h))]
    #
    # dz = SA[sol[1],sol[2]]
    nz = length(z)
    dz = sol[SVector{nz}(1:nz)]
    ∇qℓ = dz
    #
    ∇Qℓ = .5*(dz*z' + z*dz')
    #
    dℓ_dq = ∇qℓ'
    dℓ_dQv = vec(∇Qℓ)'
    #
    return dℓ_dq, dℓ_dQv
end


function proximity_jacobians(capsule_1::Capsule{T},capsule_2::Capsule{T};pdip_tol = 1e-6) where {T}

    # get the two points for each capsule
    a,b = get_ends(capsule_1)
    c,d = get_ends(capsule_2)

    # qp cost terms
    Q, q, G, h = get_cost_terms(a,b,c,d)

    # z, λ =  active_set_qp(Q,q;get_duals = true)
    if check_parallel(a,b,c,d)
        z,_,λ = pdip(Q,q,G,h; tol = pdip_tol)
    else
        z,λ = active_set_qp(Q,q;get_duals = true)
    end

    dℓ_da, dℓ_db, dℓ_dc, dℓ_dd, ∇ℓ_z = dℓ_dpoints(z,a,b,c,d)

    # optnet style QP solve + IFT for dℓ_(Qv,q)
    dℓ_dq, dℓ_dQv = optnet(z,λ,∇ℓ_z,Q,G,h)

    da_dp1, db_dp1, dc_dp2, dd_dp2 = dpoints_d_r_p(capsule_1,capsule_2)

    dQv_dr1, dQv_dp1, dQv_dr2, dQv_dp2, dq_dr1, dq_dp1, dq_dr2, dq_dp2 = dqpcost_dx(a,b,c,d,da_dp1, db_dp1, dc_dp2, dd_dp2)

    dℓ_dr1 = dℓ_da        + dℓ_db        + dℓ_dQv*dQv_dr1 + dℓ_dq*dq_dr1
    dℓ_dp1 = dℓ_da*da_dp1 + dℓ_db*db_dp1 + dℓ_dQv*dQv_dp1 + dℓ_dq*dq_dp1
    dℓ_dr2 = dℓ_dc        + dℓ_dd        + dℓ_dQv*dQv_dr2 + dℓ_dq*dq_dr2
    dℓ_dp2 = dℓ_dc*dc_dp2 + dℓ_dd*dd_dp2 + dℓ_dQv*dQv_dp2 + dℓ_dq*dq_dp2

    return dℓ_dr1, dℓ_dp1, dℓ_dr2, dℓ_dp2

end
