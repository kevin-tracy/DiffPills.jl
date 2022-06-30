function active_set_qp(Q::SMatrix{2, 2, T, 4},
                       q::SVector{2, T};
                       get_duals::Bool = false) where {T}

    # solve for unconstrained solution
    z = -Q\q

    # check if the unconstrained solution is feasible
    if in_0_1(z)
        if get_duals
            # duals are 0 if unconstrained solution is feasible
            return z, SA[0,0,0,0.0]
        else
            return z
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
            return Z[idx_min], λ
        else
            return Z[idx_min]
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
@inline function recover_duals(z::SVector{2,T},
                               Q::SMatrix{2,2,T,4},
                               q::SVector{2,T};
                               tol::T = 1e-12) where {T}

    # stationarity says Qx + q + G'z = 0, so G'z = -Qx - q
    # we are going to call y = -Qz - q to ease notation
    y = -Q*z - q
    y1,y2 = y

    # # since G' = [I(2) -I(2)], and z ≥ 0, we can recover z
    SA[
        ((y1 >= tol)  ?  y1 : 0.0 ),
        ((y2 >= tol)  ?  y2 : 0.0 ),
        ((y1 <= -tol) ? -y1 : 0.0 ),
        ((y2 <= -tol) ? -y2 : 0.0 )
    ]
end
