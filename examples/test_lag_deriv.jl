cd("/Users/kevintracy/.julia/dev/DiffPills")
import Pkg; Pkg.activate(".")
import DiffPills as dp

cd("/Users/kevintracy/.julia/dev/DiffPills/examples")
import Pkg; Pkg.activate(".")
using LinearAlgebra
using StaticArrays
using BenchmarkTools
import ForwardDiff as fd
import FiniteDiff as fd2


mutable struct CapsuleMRP{T}
    r::SVector{3,T}
    p::SVector{3,T}
    R::T
    L::T
    function CapsuleMRP(R::T,L::T) where {T}
        new{T}(SA[0,0,0.0], SA[0,0,0.0], R, L)
    end
end

@inline function get_ends(capsule::CapsuleMRP{T}) where {T}
    n_Q_b = dp.dcm_from_p(capsule.p)
    a = capsule.r + n_Q_b*SA[-capsule.L/2,0,0.0]
    b = capsule.r + n_Q_b*SA[capsule.L/2,0,0.0]
    return a, b
end
@inline function get_ends(capsule::CapsuleMRP{T},r::SVector{3,T2},p::SVector{3,T3}) where {T,T2,T3}
    n_Q_b = dp.dcm_from_p(p)
    a = r + n_Q_b*SA[-capsule.L/2,0,0.0]
    b = r + n_Q_b*SA[capsule.L/2,0,0.0]
    return a, b
end
# @inline function ℓ_fx(a::SVector{3,T},
#                       b::SVector{3,T},
#                       c::SVector{3,T},
#                       d::SVector{3,T},
#                       z::SVector{2,T},
#                       R1::T,
#                       R2::T) where {T}
#
#     z1, z2 = z
#
#     p1 = z1*a + (1-z1)*b
#     p2 = z2*c + (1-z2)*d
#
#     e = p1 - p2
#
#     # return SA[(R1 + R2)^2 - dot(e,e)]
#     return dot(e,e) - (R1 + R2)^2
# end

@inline function get_cost_terms(a::SVector{3,T},
                                b::SVector{3,T},
                                c::SVector{3,T},
                                d::SVector{3,T}) where {T}

    F = [(a-b) (d-c)]
    g = b - d

    Q = 2*F'*F
    q = 2*F'*g

    return Q, q, dot(g,g)
end

# @inline function check_parallel(a::SVector{3,T},
#                                 b::SVector{3,T},
#                                 c::SVector{3,T},
#                                 d::SVector{3,T}) where {T}
#     e1 = a - b
#     e2 = c - d
#     norm(cross(e1,e2)) < 1e-4
# end
@inline function active_set_qp(Q::SMatrix{2, 2, T, 4},
                               q::SVector{2, T};
                               get_duals::Bool = false) where {T}

    # solve for unconstrained solution
    z = -Q\q

    # check if the unconstrained solution is feasible
    if in_0_1(z)
        if get_duals
            # duals are 0 if unconstrained solution is feasible
            return z, cost(z,Q,q;check_feasibility=false), SA[0,0,0,0.0]
        else
            return z, cost(z,Q,q;check_feasibility=false)
        end
    else
        # if the unconstrained solution is infeasible, it means one
        # of the following 8 points are the solution

        # pull data from the Q and q
        Q1 = Q[1,1]
        Q2 = Q[1,2]
        Q3 = Q[2,2]
        q1,q2 = q

        # these four points are for 1 active constraint
        z1 = SA[1, -(Q2 + q2)/Q3]
        J1 = cost(z1,Q,q)

        z2 = SA[0, -q2/Q3]
        J2 = cost(z2,Q,q)

        z3 = SA[-(Q2 + q1)/Q1, 1]
        J3 = cost(z3,Q,q)

        z4 = SA[-q1/Q1, 0]
        J4 = cost(z4,Q,q)

        # these four points are if both constraints are active
        z5 = SA[0,0.0]
        J5 = 0.0 # 0 since z5 is all zeros

        z6 = SA[0,1.0]
        J6 = cost(z6,Q,q;check_feasibility=false)

        z7 = SA[1,0.0]
        J7 = cost(z7,Q,q;check_feasibility=false)

        z8 = SA[1,1.0]
        J8 = cost(z8,Q,q;check_feasibility=false)

        # find the min of these
        Z = SA[z1,z2,z3,z4,z5,z6,z7,z8]
        J = SA[J1,J2,J3,J4,J5,J6,J7,J8]

        idx_min = argmin(J)

        z = Z[idx_min]

        if get_duals
            λ = recover_duals(z,Q,q)
            return Z[idx_min], J[idx_min], λ
        else
            return Z[idx_min], J[idx_min]
        end
    end
end
@inline function in_0_1(z::SVector{2,T}) where{T}
    z1,z2 = z
    if z1<0 || z1>1
        return false
    end
    if z2<0 || z2>1
        return false
    end
    return true
end
@inline function cost(z::SVector{2,T},
                      Q::SMatrix{2, 2, T, 4},
                      q::SVector{2, T};
                      check_feasibility::Bool = true) where {T}

    if check_feasibility
        if !in_0_1(z)
            return Inf
        else
            return 0.5*z'*Q*z + q'z
        end
    else
        return 0.5*z'*Q*z + q'z
    end
end
# @inline function recover_duals(z::SVector{2,T},
#                                Q::SMatrix{2,2,T,4},
#                                q::SVector{2,T};
#                                tol::T = 1e-12) where {T}
#
#     # stationarity says Qx + q + G'z = 0, so G'z = -Qx - q
#     # we are going to call y = -Qz - q to ease notation
#     y = -Q*z - q
#     y1,y2 = y
#
#     # # since G' = [I(2) -I(2)], and z ≥ 0, we can recover z
#     SA[
#         ((y1 >= tol)  ?  y1 : 0.0 ),
#         ((y2 >= tol)  ?  y2 : 0.0 ),
#         ((y1 <= -tol) ? -y1 : 0.0 ),
#         ((y2 <= -tol) ? -y2 : 0.0 )
#     ]
# end

@inline function proximity2(P1::CapsuleMRP{T},
                           P2::CapsuleMRP{T};
                           pdip_tol = 1e-6) where {T}

    a,b = get_ends(P1)
    c,d = get_ends(P2)

    Q, q, r = get_cost_terms(a,b,c,d)

    z, J = active_set_qp(Q,q;get_duals = false)

    return J + r - (P1.R + P2.R)^2
end

# function second_qp_solver(a,b,c,d)
#     r1 = norm(a - c)
#     r2 = norm(a - d)
#     r3 = norm(b - c)
#     r4 = norm(b - d)
#
#     R = SA[r1,r2,r3,r4]
#
#     idx = argmin(R)
#
#     if idx == 1
#         return SA[0,0]
#     elseif idx == 2
#         return SA[0,]
function proximity_fd(P1,P2,r1,p1,r2,p2)
    P1.r = r1
    P1.p = p1
    P2.r = r2
    P2.p = p2
    [proximity2(P1,P2)]
end
@inline function lagrangian(P1::CapsuleMRP{T},
                            P2::CapsuleMRP{T},
                            r1::SVector{3,T2},
                            p1::SVector{3,T3},
                            r2::SVector{3,T4},
                            p2::SVector{3,T5},
                            z::SVector{2,T}) where {T,T2,T3,T4,T5}

    a,b = get_ends(P1,r1,p1)
    c,d = get_ends(P2,r2,p2)

    Q, q, r = get_cost_terms(a,b,c,d)

    0.5*z'*Q*z + q'*z + r#+ λ'*(G*z - h)
end

# @inline function lagrangian_grad(P1::CapsuleMRP{T},
#                             P2::CapsuleMRP{T},
#                             r1::SVector{3,T2},
#                             p1::SVector{3,T3},
#                             r2::SVector{3,T4},
#                             p2::SVector{3,T5},
#                             z::SVector{2,T}) where {T,T2,T3,T4,T5}
#
#     r1x,r1y,r1z = r1
#     p1x,p1y,p1z = p1
#     r2x,r2y,r2z = r2
#     p2x,p2y,p2z = p2
#
#     z1,z2 = z
#
#     R1 = P1.R
#     L1 = P1.L
#     R2 = P2.R
#     L2 = P2.L
#
# #     SA[
#     L1*z1*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1) - L2*z2*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1)
# (L1*z1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - (L2*z2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(p2x^2 + p2y^2 + p2z^2 + 1)^2
# (L2*z2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - (L1*z1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(p1x^2 + p1y^2 + p1z^2 + 1)^2
# z1*z2*((L1*L2*(8*p1y - 8*p1x*p1z)*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) - (L1*L2*(8*p1z + 8*p1x*p1y)*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (4*L1*L2*p1x*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1)*(8*p1y^2 + 8*p1z^2))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 + (4*L1*L2*p1x*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (4*L1*L2*p1x*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2)) - z1*((L1*(8*p1z + 8*p1x*p1y)*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2)))/(p1x^2 + p1y^2 + p1z^2 + 1)^2 + (L1*(8*p1y - 8*p1x*p1z)*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2)))/(p1x^2 + p1y^2 + p1z^2 + 1)^2 + (4*L1^2*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(- p1x^3*p1y - 3*p1x^2*p1z - p1x*p1y^3 - p1x*p1y*p1z^2 + 3*p1x*p1y + p1y^2*p1z + p1z^3 + p1z))/(p1x^2 + p1y^2 + p1z^2 + 1)^5 - (4*L1^2*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(p1x^3*p1z - 3*p1x^2*p1y + p1x*p1y^2*p1z + p1x*p1z^3 - 3*p1x*p1z + p1y^3 + p1y*p1z^2 + p1y))/(p1x^2 + p1y^2 + p1z^2 + 1)^5 - (2*L1^2*p1x*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1)*(8*p1y^2 + 8*p1z^2))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 - (4*L1*p1x*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 + (4*L1*p1x*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 + (4*L1*p1x*(8*p1y^2 + 8*p1z^2)*(r1x - r2x - (L1*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1))/2 + (L2*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1))/2))/(p1x^2 + p1y^2 + p1z^2 + 1)^3) - z2*((4*L1*L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z)*(p1x^3*p1z - 3*p1x^2*p1y + p1x*p1y^2*p1z + p1x*p1z^3 - 3*p1x*p1z + p1y^3 + p1y*p1z^2 + p1y))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2) - (4*L1*L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y)*(- p1x^3*p1y - 3*p1x^2*p1z - p1x*p1y^3 - p1x*p1y*p1z^2 + 3*p1x*p1y + p1y^2*p1z + p1z^3 + p1z))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (2*L1*L2*p1x*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1)*(8*p1y^2 + 8*p1z^2))/(p1x^2 + p1y^2 + p1z^2 + 1)^3)
# z2*((2*L1*L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y)*(p1x^4 + 2*p1x^2*p1z^2 - 8*p1x*p1y*p1z - p1y^4 + 6*p1y^2 + p1z^4 - 1))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2) - (4*L1*L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z)*(p1x^3 + p1x^2*p1y*p1z - 3*p1x*p1y^2 + p1x*p1z^2 + p1x + p1y^3*p1z + p1y*p1z^3 - 3*p1y*p1z))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (8*L1*L2*p1y*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/(p1x^2 + p1y^2 + p1z^2 + 1)^3) - z1*((2*L1^2*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(p1x^4 + 2*p1x^2*p1z^2 - 8*p1x*p1y*p1z - p1y^4 + 6*p1y^2 + p1z^4 - 1))/(p1x^2 + p1y^2 + p1z^2 + 1)^5 + (L1*(8*p1x - 8*p1y*p1z)*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2)))/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - (4*L1^2*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(p1x^3 + p1x^2*p1y*p1z - 3*p1x*p1y^2 + p1x*p1z^2 + p1x + p1y^3*p1z + p1y*p1z^3 - 3*p1y*p1z))/(p1x^2 + p1y^2 + p1z^2 + 1)^5 + (L1*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p1x^2 + 12*p1y^2 + 4*p1z^2 - 4))/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - (16*L1*p1y*(r1x - r2x - (L1*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1))/2 + (L2*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1))/2)*(p1x^2 - p1y^2 - p1z^2 + 1))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 - (4*L1*p1y*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 + (4*L1*p1y*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 + (8*L1^2*p1y*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/(p1x^2 + p1y^2 + p1z^2 + 1)^3) + z1*z2*((L1*L2*(8*p1x - 8*p1y*p1z)*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) - (L1*L2*(4*p1x^2 + 12*p1y^2 + 4*p1z^2 - 4)*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) - (16*L1*L2*p1y*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 + (4*L1*L2*p1y*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (4*L1*L2*p1y*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2))
# z2*((2*L1*L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z)*(p1x^4 + 2*p1x^2*p1y^2 + 8*p1x*p1y*p1z + p1y^4 - p1z^4 + 6*p1z^2 - 1))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (4*L1*L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y)*(p1x^3 - p1x^2*p1y*p1z + p1x*p1y^2 - 3*p1x*p1z^2 + p1x - p1y^3*p1z - p1y*p1z^3 + 3*p1y*p1z))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (8*L1*L2*p1z*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/(p1x^2 + p1y^2 + p1z^2 + 1)^3) - z1*((2*L1^2*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(p1x^4 + 2*p1x^2*p1y^2 + 8*p1x*p1y*p1z + p1y^4 - p1z^4 + 6*p1z^2 - 1))/(p1x^2 + p1y^2 + p1z^2 + 1)^5 + (L1*(8*p1x + 8*p1y*p1z)*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2)))/(p1x^2 + p1y^2 + p1z^2 + 1)^2 + (4*L1^2*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(p1x^3 - p1x^2*p1y*p1z + p1x*p1y^2 - 3*p1x*p1z^2 + p1x - p1y^3*p1z - p1y*p1z^3 + 3*p1y*p1z))/(p1x^2 + p1y^2 + p1z^2 + 1)^5 - (L1*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p1x^2 + 4*p1y^2 + 12*p1z^2 - 4))/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - (16*L1*p1z*(r1x - r2x - (L1*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1))/2 + (L2*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1))/2)*(p1x^2 - p1y^2 - p1z^2 + 1))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 - (4*L1*p1z*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 + (4*L1*p1z*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 + (8*L1^2*p1z*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/(p1x^2 + p1y^2 + p1z^2 + 1)^3) - z1*z2*((L1*L2*(4*p1x^2 + 4*p1y^2 + 12*p1z^2 - 4)*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (L1*L2*(8*p1x + 8*p1y*p1z)*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (16*L1*L2*p1z*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/(p1x^2 + p1y^2 + p1z^2 + 1)^3 - (4*L1*L2*p1z*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2) - (4*L1*L2*p1z*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^3*(p2x^2 + p2y^2 + p2z^2 + 1)^2))
#   L2*z2*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1) - L1*z1*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1)
# (L2*z2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - (L1*z1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(p1x^2 + p1y^2 + p1z^2 + 1)^2
# (L1*z1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - (L2*z2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(p2x^2 + p2y^2 + p2z^2 + 1)^2
# z2*((L2*(8*p2z + 8*p2x*p2y)*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2)))/(p2x^2 + p2y^2 + p2z^2 + 1)^2 + (L2*(8*p2y - 8*p2x*p2z)*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2)))/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - (4*L2^2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y)*(- p2x^3*p2y - 3*p2x^2*p2z - p2x*p2y^3 - p2x*p2y*p2z^2 + 3*p2x*p2y + p2y^2*p2z + p2z^3 + p2z))/(p2x^2 + p2y^2 + p2z^2 + 1)^5 + (4*L2^2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z)*(p2x^3*p2z - 3*p2x^2*p2y + p2x*p2y^2*p2z + p2x*p2z^3 - 3*p2x*p2z + p2y^3 + p2y*p2z^2 + p2y))/(p2x^2 + p2y^2 + p2z^2 + 1)^5 + (2*L2^2*p2x*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1)*(8*p2y^2 + 8*p2z^2))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 - (4*L2*p2x*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 + (4*L2*p2x*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 + (4*L2*p2x*(8*p2y^2 + 8*p2z^2)*(r1x - r2x - (L1*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1))/2 + (L2*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1))/2))/(p2x^2 + p2y^2 + p2z^2 + 1)^3) - z1*((4*L1*L2*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(p2x^3*p2z - 3*p2x^2*p2y + p2x*p2y^2*p2z + p2x*p2z^3 - 3*p2x*p2z + p2y^3 + p2y*p2z^2 + p2y))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3) - (4*L1*L2*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(- p2x^3*p2y - 3*p2x^2*p2z - p2x*p2y^3 - p2x*p2y*p2z^2 + 3*p2x*p2y + p2y^2*p2z + p2z^3 + p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3) + (2*L1*L2*p2x*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1)*(8*p2y^2 + 8*p2z^2))/(p2x^2 + p2y^2 + p2z^2 + 1)^3) + z1*z2*((L1*L2*(8*p2y - 8*p2x*p2z)*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) - (L1*L2*(8*p2z + 8*p2x*p2y)*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (4*L1*L2*p2x*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1)*(8*p2y^2 + 8*p2z^2))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 + (4*L1*L2*p2x*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3) + (4*L1*L2*p2x*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3))
# z1*((2*L1*L2*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(p2x^4 + 2*p2x^2*p2z^2 - 8*p2x*p2y*p2z - p2y^4 + 6*p2y^2 + p2z^4 - 1))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3) - (4*L1*L2*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(p2x^3 + p2x^2*p2y*p2z - 3*p2x*p2y^2 + p2x*p2z^2 + p2x + p2y^3*p2z + p2y*p2z^3 - 3*p2y*p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3) + (8*L1*L2*p2y*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1)*(p2x^2 - p2y^2 - p2z^2 + 1))/(p2x^2 + p2y^2 + p2z^2 + 1)^3) - z2*((2*L2^2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y)*(p2x^4 + 2*p2x^2*p2z^2 - 8*p2x*p2y*p2z - p2y^4 + 6*p2y^2 + p2z^4 - 1))/(p2x^2 + p2y^2 + p2z^2 + 1)^5 - (L2*(8*p2x - 8*p2y*p2z)*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2)))/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - (4*L2^2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z)*(p2x^3 + p2x^2*p2y*p2z - 3*p2x*p2y^2 + p2x*p2z^2 + p2x + p2y^3*p2z + p2y*p2z^3 - 3*p2y*p2z))/(p2x^2 + p2y^2 + p2z^2 + 1)^5 - (L2*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p2x^2 + 12*p2y^2 + 4*p2z^2 - 4))/(p2x^2 + p2y^2 + p2z^2 + 1)^2 + (16*L2*p2y*(r1x - r2x - (L1*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1))/2 + (L2*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1))/2)*(p2x^2 - p2y^2 - p2z^2 + 1))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 + (4*L2*p2y*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 - (4*L2*p2y*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 + (8*L2^2*p2y*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1)*(p2x^2 - p2y^2 - p2z^2 + 1))/(p2x^2 + p2y^2 + p2z^2 + 1)^3) + z1*z2*((L1*L2*(8*p2x - 8*p2y*p2z)*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) - (L1*L2*(4*p2x^2 + 12*p2y^2 + 4*p2z^2 - 4)*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) - (16*L1*L2*p2y*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1)*(p2x^2 - p2y^2 - p2z^2 + 1))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 + (4*L1*L2*p2y*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3) + (4*L1*L2*p2y*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3))
# z1*((2*L1*L2*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(p2x^4 + 2*p2x^2*p2y^2 + 8*p2x*p2y*p2z + p2y^4 - p2z^4 + 6*p2z^2 - 1))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3) + (4*L1*L2*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(p2x^3 - p2x^2*p2y*p2z + p2x*p2y^2 - 3*p2x*p2z^2 + p2x - p2y^3*p2z - p2y*p2z^3 + 3*p2y*p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3) + (8*L1*L2*p2z*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1)*(p2x^2 - p2y^2 - p2z^2 + 1))/(p2x^2 + p2y^2 + p2z^2 + 1)^3) - z2*((2*L2^2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z)*(p2x^4 + 2*p2x^2*p2y^2 + 8*p2x*p2y*p2z + p2y^4 - p2z^4 + 6*p2z^2 - 1))/(p2x^2 + p2y^2 + p2z^2 + 1)^5 - (L2*(8*p2x + 8*p2y*p2z)*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2)))/(p2x^2 + p2y^2 + p2z^2 + 1)^2 + (4*L2^2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y)*(p2x^3 - p2x^2*p2y*p2z + p2x*p2y^2 - 3*p2x*p2z^2 + p2x - p2y^3*p2z - p2y*p2z^3 + 3*p2y*p2z))/(p2x^2 + p2y^2 + p2z^2 + 1)^5 + (L2*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p2x^2 + 4*p2y^2 + 12*p2z^2 - 4))/(p2x^2 + p2y^2 + p2z^2 + 1)^2 + (16*L2*p2z*(r1x - r2x - (L1*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1))/2 + (L2*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1))/2)*(p2x^2 - p2y^2 - p2z^2 + 1))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 + (4*L2*p2z*(r1z - r2z + (L1*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) - (L2*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 - (4*L2*p2z*(r1y - r2y - (L1*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/(2*(p1x^2 + p1y^2 + p1z^2 + 1)^2) + (L2*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(2*(p2x^2 + p2y^2 + p2z^2 + 1)^2))*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 + (8*L2^2*p2z*((8*p2y^2 + 8*p2z^2)/(p2x^2 + p2y^2 + p2z^2 + 1)^2 - 1)*(p2x^2 - p2y^2 - p2z^2 + 1))/(p2x^2 + p2y^2 + p2z^2 + 1)^3) - z1*z2*((L1*L2*(4*p2x^2 + 4*p2y^2 + 12*p2z^2 - 4)*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (L1*L2*(8*p2x + 8*p2y*p2z)*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^2) + (16*L1*L2*p2z*((8*p1y^2 + 8*p1z^2)/(p1x^2 + p1y^2 + p1z^2 + 1)^2 - 1)*(p2x^2 - p2y^2 - p2z^2 + 1))/(p2x^2 + p2y^2 + p2z^2 + 1)^3 - (4*L1*L2*p2z*(4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)*(4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3) - (4*L1*L2*p2z*(4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)*(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z))/((p1x^2 + p1y^2 + p1z^2 + 1)^2*(p2x^2 + p2y^2 + p2z^2 + 1)^3))
# ]
#
# first_term = (p2x^2 + p2y^2 + p2z^2 + 1)
# second_term = (p1x^2 + p1y^2 + p1z^2 + 1)
# third_term=(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z)
# fourth_term = (4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)
# fifth_term = (4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y)
# sixth_term = (8*p2y^2 + 8*p2z^2)
# seventh_term = (8*p1y^2 + 8*p1z^2)
# eighth_term = (p1x^4 + 2*p1x^2*p1z^2 - 8*p1x*p1y*p1z - p1y^4 + 6*p1y^2 + p1z^4 - 1)
# ninth_term = (4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)
# tenth_term = (p2x^2 - p2y^2 - p2z^2 + 1)
# eleventh_term = (p2x^3 - p2x^2*p2y*p2z + p2x*p2y^2 - 3*p2x*p2z^2 + p2x - p2y^3*p2z - p2y*p2z^3 + 3*p2y*p2z)
#
# SA[
# L1*z1*(seventh_term/second_term^2 - 1) - L2*z2*(sixth_term/first_term^2 - 1)
# (L1*z1*ninth_term)/second_term^2 - (L2*z2*third_term)/first_term^2
# (L2*z2*fifth_term)/first_term^2 - (L1*z1*fourth_term)/second_term^2
# z1*z2*((L1*L2*(8*p1y - 8*p1x*p1z)*third_term)/(second_term^2*first_term^2) - (L1*L2*(8*p1z + 8*p1x*p1y)*fifth_term)/(second_term^2*first_term^2) + (4*L1*L2*p1x*(sixth_term/first_term^2 - 1)*seventh_term)/second_term^3 + (4*L1*L2*p1x*fourth_term*fifth_term)/(second_term^3*first_term^2) + (4*L1*L2*p1x*ninth_term*third_term)/(second_term^3*first_term^2)) - z1*((L1*(8*p1z + 8*p1x*p1y)*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2)))/second_term^2 + (L1*(8*p1y - 8*p1x*p1z)*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2)))/second_term^2 + (4*L1^2*fourth_term*(- p1x^3*p1y - 3*p1x^2*p1z - p1x*p1y^3 - p1x*p1y*p1z^2 + 3*p1x*p1y + p1y^2*p1z + p1z^3 + p1z))/second_term^5 - (4*L1^2*ninth_term*(p1x^3*p1z - 3*p1x^2*p1y + p1x*p1y^2*p1z + p1x*p1z^3 - 3*p1x*p1z + p1y^3 + p1y*p1z^2 + p1y))/second_term^5 - (2*L1^2*p1x*(seventh_term/second_term^2 - 1)*seventh_term)/second_term^3 - (4*L1*p1x*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fourth_term)/second_term^3 + (4*L1*p1x*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*ninth_term)/second_term^3 + (4*L1*p1x*seventh_term*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2))/second_term^3) - z2*((4*L1*L2*third_term*(p1x^3*p1z - 3*p1x^2*p1y + p1x*p1y^2*p1z + p1x*p1z^3 - 3*p1x*p1z + p1y^3 + p1y*p1z^2 + p1y))/(second_term^3*first_term^2) - (4*L1*L2*fifth_term*(- p1x^3*p1y - 3*p1x^2*p1z - p1x*p1y^3 - p1x*p1y*p1z^2 + 3*p1x*p1y + p1y^2*p1z + p1z^3 + p1z))/(second_term^3*first_term^2) + (2*L1*L2*p1x*(sixth_term/first_term^2 - 1)*seventh_term)/second_term^3)
# z2*((2*L1*L2*fifth_term*eighth_term)/(second_term^3*first_term^2) - (4*L1*L2*third_term*(p1x^3 + p1x^2*p1y*p1z - 3*p1x*p1y^2 + p1x*p1z^2 + p1x + p1y^3*p1z + p1y*p1z^3 - 3*p1y*p1z))/(second_term^3*first_term^2) + (8*L1*L2*p1y*(sixth_term/first_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3) - z1*((2*L1^2*fourth_term*eighth_term)/second_term^5 + (L1*(8*p1x - 8*p1y*p1z)*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2)))/second_term^2 - (4*L1^2*ninth_term*(p1x^3 + p1x^2*p1y*p1z - 3*p1x*p1y^2 + p1x*p1z^2 + p1x + p1y^3*p1z + p1y*p1z^3 - 3*p1y*p1z))/second_term^5 + (L1*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*(4*p1x^2 + 12*p1y^2 + 4*p1z^2 - 4))/second_term^2 - (16*L1*p1y*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3 - (4*L1*p1y*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fourth_term)/second_term^3 + (4*L1*p1y*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*ninth_term)/second_term^3 + (8*L1^2*p1y*(seventh_term/second_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3) + z1*z2*((L1*L2*(8*p1x - 8*p1y*p1z)*third_term)/(second_term^2*first_term^2) - (L1*L2*(4*p1x^2 + 12*p1y^2 + 4*p1z^2 - 4)*fifth_term)/(second_term^2*first_term^2) - (16*L1*L2*p1y*(sixth_term/first_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3 + (4*L1*L2*p1y*fourth_term*fifth_term)/(second_term^3*first_term^2) + (4*L1*L2*p1y*ninth_term*third_term)/(second_term^3*first_term^2))
# z2*((2*L1*L2*third_term*(p1x^4 + 2*p1x^2*p1y^2 + 8*p1x*p1y*p1z + p1y^4 - p1z^4 + 6*p1z^2 - 1))/(second_term^3*first_term^2) + (4*L1*L2*fifth_term*(p1x^3 - p1x^2*p1y*p1z + p1x*p1y^2 - 3*p1x*p1z^2 + p1x - p1y^3*p1z - p1y*p1z^3 + 3*p1y*p1z))/(second_term^3*first_term^2) + (8*L1*L2*p1z*(sixth_term/first_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3) - z1*((2*L1^2*ninth_term*(p1x^4 + 2*p1x^2*p1y^2 + 8*p1x*p1y*p1z + p1y^4 - p1z^4 + 6*p1z^2 - 1))/second_term^5 + (L1*(8*p1x + 8*p1y*p1z)*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2)))/second_term^2 + (4*L1^2*fourth_term*(p1x^3 - p1x^2*p1y*p1z + p1x*p1y^2 - 3*p1x*p1z^2 + p1x - p1y^3*p1z - p1y*p1z^3 + 3*p1y*p1z))/second_term^5 - (L1*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*(4*p1x^2 + 4*p1y^2 + 12*p1z^2 - 4))/second_term^2 - (16*L1*p1z*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3 - (4*L1*p1z*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fourth_term)/second_term^3 + (4*L1*p1z*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*ninth_term)/second_term^3 + (8*L1^2*p1z*(seventh_term/second_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3) - z1*z2*((L1*L2*(4*p1x^2 + 4*p1y^2 + 12*p1z^2 - 4)*third_term)/(second_term^2*first_term^2) + (L1*L2*(8*p1x + 8*p1y*p1z)*fifth_term)/(second_term^2*first_term^2) + (16*L1*L2*p1z*(sixth_term/first_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3 - (4*L1*L2*p1z*fourth_term*fifth_term)/(second_term^3*first_term^2) - (4*L1*L2*p1z*ninth_term*third_term)/(second_term^3*first_term^2))
# L2*z2*(sixth_term/first_term^2 - 1) - L1*z1*(seventh_term/second_term^2 - 1)
# (L2*z2*third_term)/first_term^2 - (L1*z1*ninth_term)/second_term^2
# (L1*z1*fourth_term)/second_term^2 - (L2*z2*fifth_term)/first_term^2
# z2*((L2*(8*p2z + 8*p2x*p2y)*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2)))/first_term^2 + (L2*(8*p2y - 8*p2x*p2z)*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2)))/first_term^2 - (4*L2^2*fifth_term*(- p2x^3*p2y - 3*p2x^2*p2z - p2x*p2y^3 - p2x*p2y*p2z^2 + 3*p2x*p2y + p2y^2*p2z + p2z^3 + p2z))/first_term^5 + (4*L2^2*third_term*(p2x^3*p2z - 3*p2x^2*p2y + p2x*p2y^2*p2z + p2x*p2z^3 - 3*p2x*p2z + p2y^3 + p2y*p2z^2 + p2y))/first_term^5 + (2*L2^2*p2x*(sixth_term/first_term^2 - 1)*sixth_term)/first_term^3 - (4*L2*p2x*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fifth_term)/first_term^3 + (4*L2*p2x*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*third_term)/first_term^3 + (4*L2*p2x*sixth_term*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2))/first_term^3) - z1*((4*L1*L2*ninth_term*(p2x^3*p2z - 3*p2x^2*p2y + p2x*p2y^2*p2z + p2x*p2z^3 - 3*p2x*p2z + p2y^3 + p2y*p2z^2 + p2y))/(second_term^2*first_term^3) - (4*L1*L2*fourth_term*(- p2x^3*p2y - 3*p2x^2*p2z - p2x*p2y^3 - p2x*p2y*p2z^2 + 3*p2x*p2y + p2y^2*p2z + p2z^3 + p2z))/(second_term^2*first_term^3) + (2*L1*L2*p2x*(seventh_term/second_term^2 - 1)*sixth_term)/first_term^3) + z1*z2*((L1*L2*(8*p2y - 8*p2x*p2z)*ninth_term)/(second_term^2*first_term^2) - (L1*L2*(8*p2z + 8*p2x*p2y)*fourth_term)/(second_term^2*first_term^2) + (4*L1*L2*p2x*(seventh_term/second_term^2 - 1)*sixth_term)/first_term^3 + (4*L1*L2*p2x*fourth_term*fifth_term)/(second_term^2*first_term^3) + (4*L1*L2*p2x*ninth_term*third_term)/(second_term^2*first_term^3))
# z1*((2*L1*L2*fourth_term*(p2x^4 + 2*p2x^2*p2z^2 - 8*p2x*p2y*p2z - p2y^4 + 6*p2y^2 + p2z^4 - 1))/(second_term^2*first_term^3) - (4*L1*L2*ninth_term*(p2x^3 + p2x^2*p2y*p2z - 3*p2x*p2y^2 + p2x*p2z^2 + p2x + p2y^3*p2z + p2y*p2z^3 - 3*p2y*p2z))/(second_term^2*first_term^3) + (8*L1*L2*p2y*(seventh_term/second_term^2 - 1)*tenth_term)/first_term^3) - z2*((2*L2^2*fifth_term*(p2x^4 + 2*p2x^2*p2z^2 - 8*p2x*p2y*p2z - p2y^4 + 6*p2y^2 + p2z^4 - 1))/first_term^5 - (L2*(8*p2x - 8*p2y*p2z)*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2)))/first_term^2 - (4*L2^2*third_term*(p2x^3 + p2x^2*p2y*p2z - 3*p2x*p2y^2 + p2x*p2z^2 + p2x + p2y^3*p2z + p2y*p2z^3 - 3*p2y*p2z))/first_term^5 - (L2*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*(4*p2x^2 + 12*p2y^2 + 4*p2z^2 - 4))/first_term^2 + (16*L2*p2y*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2)*tenth_term)/first_term^3 + (4*L2*p2y*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fifth_term)/first_term^3 - (4*L2*p2y*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*third_term)/first_term^3 + (8*L2^2*p2y*(sixth_term/first_term^2 - 1)*tenth_term)/first_term^3) + z1*z2*((L1*L2*(8*p2x - 8*p2y*p2z)*ninth_term)/(second_term^2*first_term^2) - (L1*L2*(4*p2x^2 + 12*p2y^2 + 4*p2z^2 - 4)*fourth_term)/(second_term^2*first_term^2) - (16*L1*L2*p2y*(seventh_term/second_term^2 - 1)*tenth_term)/first_term^3 + (4*L1*L2*p2y*fourth_term*fifth_term)/(second_term^2*first_term^3) + (4*L1*L2*p2y*ninth_term*third_term)/(second_term^2*first_term^3))
# z1*((2*L1*L2*ninth_term*(p2x^4 + 2*p2x^2*p2y^2 + 8*p2x*p2y*p2z + p2y^4 - p2z^4 + 6*p2z^2 - 1))/(second_term^2*first_term^3) + (4*L1*L2*fourth_term*eleventh_term)/(second_term^2*first_term^3) + (8*L1*L2*p2z*(seventh_term/second_term^2 - 1)*tenth_term)/first_term^3) - z2*((2*L2^2*third_term*(p2x^4 + 2*p2x^2*p2y^2 + 8*p2x*p2y*p2z + p2y^4 - p2z^4 + 6*p2z^2 - 1))/first_term^5 - (L2*(8*p2x + 8*p2y*p2z)*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2)))/first_term^2 + (4*L2^2*fifth_term*eleventh_term)/first_term^5 + (L2*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*(4*p2x^2 + 4*p2y^2 + 12*p2z^2 - 4))/first_term^2 + (16*L2*p2z*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2)*tenth_term)/first_term^3 + (4*L2*p2z*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fifth_term)/first_term^3 - (4*L2*p2z*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*third_term)/first_term^3 + (8*L2^2*p2z*(sixth_term/first_term^2 - 1)*tenth_term)/first_term^3) - z1*z2*((L1*L2*(4*p2x^2 + 4*p2y^2 + 12*p2z^2 - 4)*ninth_term)/(second_term^2*first_term^2) + (L1*L2*(8*p2x + 8*p2y*p2z)*fourth_term)/(second_term^2*first_term^2) + (16*L1*L2*p2z*(seventh_term/second_term^2 - 1)*tenth_term)/first_term^3 - (4*L1*L2*p2z*fourth_term*fifth_term)/(second_term^2*first_term^3) - (4*L1*L2*p2z*ninth_term*third_term)/(second_term^2*first_term^3))
# ]
#
# end
# let
#
#     P1 = CapsuleMRP(0.2,0.9)
#     P2 = CapsuleMRP(0.3,0.8)
#
#     P1.r = @SVector randn(3)
#     P1.p = @SVector randn(3)
#     P2.r = @SVector randn(3)
#     P2.p = @SVector randn(3)
#
#     @show proximity(P1,P2)
#
#     @btime proximity($P1,$P2)
#
#     idx_r1 = SVector{3}(1:3)
#     idx_p1 = SVector{3}(4:6)
#     idx_r2 = SVector{3}(7:9)
#     idx_p2 = SVector{3}(10:12)
#
    # J1 =  fd2.finite_difference_jacobian(_θ -> proximity_fd(P1,P2,_θ[idx_r1],_θ[idx_p1],_θ[idx_r2],_θ[idx_p2]), [P1.r;P1.p;P2.r;P2.p])
#
#     a,b = get_ends(P1)
#     c,d = get_ends(P2)
#
#     Q, q = get_cost_terms(a,b,c,d)
#
#     z = active_set_qp(Q,q;get_duals = false)
#
    # J2 = fd.gradient(_θ -> lagrangian(P1,P2,_θ[idx_r1],_θ[idx_p1],_θ[idx_r2],_θ[idx_p2],z), [P1.r;P1.p;P2.r;P2.p])
#
#     @btime fd.gradient(_θ -> lagrangian($P1,$P2,_θ[$idx_r1],_θ[$idx_p1],_θ[$idx_r2],_θ[$idx_p2],$z), [$P1.r;$P1.p;$P2.r;$P2.p])
    # @show norm(vec(J1) - vec(J2))
#
#     @btime lagrangian_grad($P1,$P2,$P1.r, $P1.p, $P2.r, $P2.p, $z )
#
#     J3 = lagrangian_grad(P1,P2,P1.r, P1.p, P2.r, P2.p, z )
#
#     @show norm(vec(J1) - vec(J3))
# end

let

    P1 = CapsuleMRP(1.4,1.8)
    P2 = CapsuleMRP(0.7,1.1)

    P1.r = SVector(2.1,-3.3,1.4)
    P1.p = SVector(0.1,0.3,0.4)
    P2.r = SVector(-2.1,-4.3,-4.4)
    P2.p = SVector(-0.23,0.11,-0.32)

    a,b = get_ends(P1)
    c,d = get_ends(P2)

    Q, q = get_cost_terms(a,b,c,d)

    z, J = active_set_qp(Q,q;get_duals = false)

        idx_r1 = SVector{3}(1:3)
        idx_p1 = SVector{3}(4:6)
        idx_r2 = SVector{3}(7:9)
        idx_p2 = SVector{3}(10:12)
    #
        J1 =  fd2.finite_difference_jacobian(_θ -> proximity_fd(P1,P2,_θ[idx_r1],_θ[idx_p1],_θ[idx_r2],_θ[idx_p2]), [P1.r;P1.p;P2.r;P2.p])
    #
    #     a,b = get_ends(P1)
    #     c,d = get_ends(P2)
    #
    #     Q, q = get_cost_terms(a,b,c,d)
    #
    #     z = active_set_qp(Q,q;get_duals = false)
    #
        J2 = fd.gradient(_θ -> lagrangian(P1,P2,_θ[idx_r1],_θ[idx_p1],_θ[idx_r2],_θ[idx_p2],z), [P1.r;P1.p;P2.r;P2.p])
    #
    #     @btime fd.gradient(_θ -> lagrangian($P1,$P2,_θ[$idx_r1],_θ[$idx_p1],_θ[$idx_r2],_θ[$idx_p2],$z), [$P1.r;$P1.p;$P2.r;$P2.p])
        @show norm(vec(J1) - vec(J2))
        @show J1
        @show J2
    #
    #     @btime lagrangian_grad($P1,$P2,$P1.r, $P1.p, $P2.r, $P2.p, $z )
    #
    #     J3 = lagrangian_grad(P1,P2,P1.r, P1.p, P2.r, P2.p, z )
    #
    #     @show norm(vec(J1) - vec(J3))

    # true_cost =  J + r
    #
    # @show true_cost
    #
    # # # @show a
    # # # @show b
    # # # @show c
    # # # @show d
    # #
    # # # Q,q,r = get_cost_terms(a,b,c,d)
    # # F = [(a-b) (d-c)]
    # # g = b - d
    # #
    # # Q = 2*F'*F
    # # q = 2*F'*g
    # # r = dot(g,g)
    # #
    # # z = randn(2)
    # p1 = z[1]*a + (1 - z[1])*b
    # p2 = z[2]*c + (1 - z[2])*d
    # e = p1 - p2
    # # # e2 = F*z + g
    # # # @show e
    # # # @show e2
    # # #
    # @show dot(e,e)
    # # # @show z'*(F'*F)*z + dot(2*F'*g,z) + dot(g,g)
    # # # @show 0.5*z'*(2*F'*F)*z + dot(2*F'*g,z) + dot(g,g)
    # # J1 = dot(e,e)
    # #
    # @show 0.5*z'*Q*z + q'*z + r
    #
    # @show proximity(P1,P2)
    # #
    # # @show abs(J1-J2)
    # # # @show dot(g,g)

end
