
mutable struct Capsule{T}
    r::SVector{3,T}
    p::SVector{3,T}
    q::SVector{4,T}
    R::T
    L::T
    attitude::Symbol
end

function create_capsule(attitude::Symbol)
    return Capsule(SA[0.0,0.0,0.0],SA[0,0,0.0],SA[1,0,0,0.0],0.0,0.0,attitude)
end

@inline function get_ends(capsule::Capsule{T}) where {T}
    if capsule.attitude == :MRP
        n_Q_b = dcm_from_p(capsule.p)
    end
    if capsule.attitude == :quat
        n_Q_b = dcm_from_q(capsule.q)
    end
    a = capsule.r + n_Q_b*SA[-capsule.L/2,0,0]
    b = capsule.r + n_Q_b*SA[capsule.L/2,0,0]
    return a, b
end

@inline function ℓ_fx(a::SVector{3,T},
                      b::SVector{3,T},
                      c::SVector{3,T},
                      d::SVector{3,T},
                      z::SVector{2,T},
                      R1::T,
                      R2::T) where {T}

    z1, z2 = z

    p1 = z1*a + (1-z1)*b
    p2 = z2*c + (1-z2)*d

    e = p1 - p2

    # return SA[(R1 + R2)^2 - dot(e,e)]
    return dot(e,e) - (R1 + R2)^2
end

@inline function get_cost_terms(a::SVector{3,T},
                                b::SVector{3,T},
                                c::SVector{3,T},
                                d::SVector{3,T}) where {T}

    F = [(a-b) (d-c)]

    Q = F'*F
    q = F'*(b-d)

    G = SA[1 0; 0 1; -1 0; 0 -1.0]
    h = SA[1,1,0,0.0]

    return Q, q, G, h
end

@inline function check_parallel(a::SVector{3,T},
                                b::SVector{3,T},
                                c::SVector{3,T},
                                d::SVector{3,T}) where {T}
    e1 = a - b
    e2 = c - d
    norm(cross(e1,e2)) < 1e-4
end

@inline function proximity(P1::Capsule{T},
                                P2::Capsule{T};
                                pdip_tol = 1e-6) where {T}

    a,b = get_ends(P1)
    c,d = get_ends(P2)

    Q, q, G, h = get_cost_terms(a,b,c,d)

    if check_parallel(a,b,c,d)
        z,_,_ = pdip(Q,q,G,h; tol = pdip_tol)
    else
        z = active_set_qp(Q,q;get_duals = false)
    end

    return ℓ_fx(a,b,c,d,z,P1.R,P2.R)
end

@inline function closest_points(P1::Capsule{T},
                                P2::Capsule{T};
                                pdip_tol = 1e-6) where {T}

    a,b = get_ends(P1)
    c,d = get_ends(P2)

    Q, q, G, h = get_cost_terms(a,b,c,d)

    if check_parallel(a,b,c,d)
        z,_,_ = pdip(Q,q,G,h; tol = pdip_tol)
    else
        z = active_set_qp(Q,q;get_duals = false)
    end

    p1 = z[1]*(a-b) + b
    p2 = z[2]*(c-d) + d

    n1 = normalize(p2 - p1)
    n2 = normalize(p1 - p2)
    # n2 = -n1 # this is also correct

    p1_edge = p1 + P1.R*n1
    p2_edge = p2 + P2.R*n2

    # edge points are always on the edge
    # n1 is the contact force direction on body 1
    return p1_edge, p2_edge, -n1, -n2
end
