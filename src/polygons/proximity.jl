
mutable struct Polygon{T,nh,nh2}
    r::SVector{3,T}         # position in world frame
    p::SVector{3,T}         # MRP for attitude (ᴺpᴮ)
    q::SVector{4,T}         # quat for attitude (ᴺqᴮ)
    A::SMatrix{nh,2,T,nh2}  # polygon description (Ax≦b)
    b::SVector{nh,T}        # polygon description (Ax≦b)
    R::T                    # "cushion" radius
    attitude::Symbol        # attitude parametrization (:MRP or :quat)
end

function create_polygon(A::SMatrix{nh,2,T,nh2},b::SVector{nh,T},R::T,attitude::Symbol) where {nh,T,nh2}
    return Polygon{T,nh,nh2}(SA[0,0,0.0],SA[0,0,0.0],SA[1,0,0,0.0],A,b,R,attitude)
end


@inline function ℓ_fx(r1::SVector{3,T},
                      r2::SVector{3,T},
                      R1::T,
                      R2::T,
                      Q̃1::SMatrix{3,2,T},
                      Q̃2::SMatrix{3,2,T},
                      z::SVector{4,T}) where {T}

    # two closest points on the two planes
    p1 = r1 + Q̃1*SA[z[1],z[2]]
    p2 = r2 + Q̃2*SA[z[3],z[4]]

    # vector between them
    e = p1 - p2

    # return the squared distance
    dot(e,e) - (R1 + R2)^2
end

@inline function get_qp_terms(P1::Polygon{T,nh_1,nh2_1},
                              P2::Polygon{T,nh_2,nh2_2},
                              Q̃1::SMatrix{3,2,T,6},
                              Q̃2::SMatrix{3,2,T,6}) where {T,nh_1,nh_2,nh2_1,nh2_2}

    F = [Q̃1 -Q̃2]

    # if the attitudes are the same, set Q manually
    if norm(Q̃1 - Q̃2) < 1e-10
        Q = SA[1  0 -1  0;
               0  1  0 -1;
              -1  0  1  0;
               0 -1  0  1.0]
    else
        Q = F'*F
    end

    q = F'*(P1.r-P2.r)

    # @show Q
    # @show eigvals(Q)

    G1 = hcat(P1.A, (@SMatrix zeros(nh_1,2)))
    G2 = hcat( (@SMatrix zeros(nh_2,2)), P2.A)
    G = vcat(G1,G2)
    h = vcat(P1.b,P2.b)
    return Q, q, G, h
end

@inline function proximity(P1::Polygon{T,nh_1,nh2_1},
                           P2::Polygon{T,nh_2,nh2_2};
                           pdip_tol = 1e-6) where {T,nh_1,nh_2,nh2_1,nh2_2}

    if P1.attitude == :MRP
        Q̃1 = dcm_tilde_from_p(P1.p)
        Q̃2 = dcm_tilde_from_p(P2.p)
    elseif P1.attitude == :quat
        Q̃1 = dcm_tilde_from_q(P1.q)
        Q̃2 = dcm_tilde_from_q(P2.q)
    end

    Q, q, G, h = get_qp_terms(P1,P2,Q̃1,Q̃2)

    z,_,_ = pdip(Q,q,G,h; tol = pdip_tol)

    return ℓ_fx(P1.r,P2.r,P1.R,P2.R,Q̃1,Q̃2,z)
end

@inline function closest_points(P1::Polygon{T,nh_1,nh2_1},
                                P2::Polygon{T,nh_2,nh2_2};
                                pdip_tol = 1e-6) where {T,nh_1,nh_2,nh2_1,nh2_2}

    if P1.attitude == :MRP
        Q̃1 = dcm_tilde_from_p(P1.p)
        Q̃2 = dcm_tilde_from_p(P2.p)
    elseif P1.attitude == :quat
        Q̃1 = dcm_tilde_from_q(P1.q)
        Q̃2 = dcm_tilde_from_q(P2.q)
    end

    Q, q, G, h = get_qp_terms(P1,P2,Q̃1,Q̃2)

    z,_,_ = pdip(Q,q,G,h; tol = pdip_tol)

    p1 = P1.r + Q̃1*SA[z[1],z[2]]
    p2 = P2.r + Q̃2*SA[z[3],z[4]]

    n1 = normalize(p2 - p1)
    n2 = normalize(p1 - p2)
    # n2 = -n1 # this is also correct

    p1_edge = p1 + P1.R*n1
    p2_edge = p2 + P2.R*n2

    # edge points are always on the edge
    # n1 is the contact force direction on body 1
    return p1_edge, p2_edge, -n1, -n2
end
