
@inline function cvp_qp_terms(a1::SVector{3,T},
                              b1::SVector{3,T},
                              Q̃::SMatrix{3,2,T,6},
                              r::SVector{3,T},
                              polygon::Polygon{T,nh,nh2}) where {T,nh,nh2}

    F = hcat((a1-b1),-Q̃)
    g = b1 - r
    Q = F'*F
    q = F'*g

    G1 = SA[ 1 0 0;
          -1 0 0.0]
    h1 = SA[1,0.0]

    G2 = hcat((@SVector zeros(nh)), polygon.A)
    h2 = polygon.b

    G = vcat(G1,G2)
    h = vcat(h1,h2)

    return Q, q, G, h
end

@inline function ℓ_cvp(z::SVector{3,T},
                       a1::SVector{3,T},
                       b1::SVector{3,T},
                       r::SVector{3,T},
                       Q̃::SMatrix{3,2,T,6},
                       cap::Capsule{T},
                       polygon::Polygon{T,nh,nh2}) where {T,nh,nh2}

    θ = z[1]
    x = z[SA[2,3]]

    p1 = θ*a1 + (1-θ)*b1
    p2 = r + Q̃*x

    e = p1 - p2

    return dot(e,e) - (cap.R + polygon.R)^2
end
@inline function proximity(cap::Capsule{T},
                           polygon::Polygon{T,nh,nh2};
                           pdip_tol::T=1e-6) where {T,nh,nh2}

    # end points of capsule
    a1,b1 = get_ends(cap)

    # polygon attitude
    if polygon.attitude == :MRP
        Q̃ = dcm_tilde_from_p(polygon.p)
    elseif polygon.attitude == :quat
        Q̃ = dcm_tilde_from_q(polygon.q)
    end

    # QP terms
    Q, q, G, h = cvp_qp_terms(a1,b1,Q̃,polygon.r, polygon)

    z,_,_ = pdip(Q,q,G,h;tol = pdip_tol)

    ℓ_cvp(z, a1, b1, polygon.r, Q̃, cap, polygon)
end

@inline function closest_points(cap::Capsule{T},
                                polygon::Polygon{T,nh,nh2};
                                pdip_tol::T=1e-6) where {T,nh,nh2}

    # end points of capsule
    a1,b1 = get_ends(cap)

    # polygon attitude
    if polygon.attitude == :MRP
        Q̃ = dcm_tilde_from_p(polygon.p)
    elseif polygon.attitude == :quat
        Q̃ = dcm_tilde_from_q(polygon.q)
    end

    # QP terms
    Q, q, G, h = cvp_qp_terms(a1,b1,Q̃,polygon.r, polygon)

    z,_,_ = pdip(Q,q,G,h;tol = pdip_tol)

    θ = z[1]
    x = z[SA[2,3]]

    p1 = θ*a1 + (1-θ)*b1
    p2 = polygon.r + Q̃*x

    n1 = normalize(p2 - p1)
    n2 = normalize(p1 - p2)

    p1_edge = p1 + cap.R*n1
    p2_edge = p2 + polygon.R*n2

    return p1_edge, p2_edge
end
# @inline function proximity(polygon::Polygon{T,nh,nh2},
#                            cap::Capsule{T};
#                            pdip_tol::T=1e-6) where {T,nh,nh2}
#     # switched order
#     proximity(cap,polygon;pdip_tol = pdip_tol)
# end
