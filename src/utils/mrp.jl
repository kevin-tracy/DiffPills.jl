@inline function dcm_from_p(p::SVector{3,T}) where {T}
    p1,p2,p3 = p
    den = (p1^2 + p2^2 + p3^2 + 1)^2
    a = (4*p1^2 + 4*p2^2 + 4*p3^2 - 4)
    Q = SA[
    (-((8*p2^2+8*p3^2)/den-1)*den)   (8*p1*p2 + p3*a)     (8*p1*p3 - p2*a);
    (8*p1*p2 - p3*a) (-((8*p1^2 + 8*p3^2)/den - 1)*den)   (8*p2*p3 + p1*a);
    (8*p1*p3 + p2*a)  (8*p2*p3 - p1*a)  (-((8*p1^2 + 8*p2^2)/den - 1)*den)
    ]/den
end

@inline function skew(p::SVector{3,T}) where {T}
    # cross product matrix
    SA[0 -p[3] p[2]; p[3] 0 -p[1]; -p[2] p[1] 0]
end

@inline function mrp_dcm_jacobian(p::SVector{3,T},v::SVector{3,T}) where {T}
    # jacobian of (dcm_from_p(p)*v) wrt p
    p1 = 8*((p'*v)*p - (p'p)*v)
    p2 = 4*(1 - p'p)*cross(p,v)
    p3 = (1 + p'p)^2
    dp1 = 8*(p*v' + (p'v)*I - 2*v*p')
    dp2 = -4*(cross(v,p)*(-2*p)' + (1 - p'p)*skew(v))
    dp3 = 4*(1 + p'p)*p'
    ((dp1 + dp2)*p3 - (p1 + p2)*dp3)/(p3^2)
end
@inline function mrp_dcm_T_jacobian(p::SVector{3,T},v::SVector{3,T}) where {T}
    # jacobian of (dcm_from_p(p)'*v) wrt p
    p = -p
    p1 = 8*((p'*v)*p - (p'p)*v)
    p2 = 4*(1 - p'p)*cross(p,v)
    p3 = (1 + p'p)^2
    dp1 = 8*(p*v' + (p'v)*I - 2*v*p')
    dp2 = -4*(cross(v,p)*(-2*p)' + (1 - p'p)*skew(v))
    dp3 = 4*(1 + p'p)*p'
    -((dp1 + dp2)*p3 - (p1 + p2)*dp3)/(p3^2)
end
@inline function mrp_dcm_tilde_jacobian(p::SVector{3,T},v::SVector{2,T}) where {T}
    mrp_dcm_jacobian(p,SA[v[1],v[2],0.0])
end
@inline function mrp_dcm_tilde_T_jacobian(p::SVector{3,T},v::SVector{3,T}) where {T}
    J = mrp_dcm_T_jacobian(p,v)
    J[SA[1,2],:]
end
@inline function dcm_tilde_from_p(p::SVector{3,T}) where {T}
    # only return first 2 columns
    p1,p2,p3 = p
    den = (p1^2 + p2^2 + p3^2 + 1)^2
    a = (4*p1^2 + 4*p2^2 + 4*p3^2 - 4)
    Q = SA[
    (-((8*p2^2+8*p3^2)/den-1)*den)   (8*p1*p2 + p3*a)     ;
    (8*p1*p2 - p3*a) (-((8*p1^2 + 8*p3^2)/den - 1)*den)   ;
    (8*p1*p3 + p2*a)  (8*p2*p3 - p1*a)
    ]/den
end

@inline function dcm_prod_jac_p1(p1::SVector{3,T},Q2::SMatrix{3,3,T,9}) where {T}
    # jacobian of vec(Q1'*Q2) wrt p1
    vcat(mrp_dcm_T_jacobian(p1,Q2[:,1]),
         mrp_dcm_T_jacobian(p1,Q2[:,2]),
         mrp_dcm_T_jacobian(p1,Q2[:,3]))
end

@inline function dcm_prod_jac_p2(Q1::SMatrix{3,3,T,9},p2::SVector{3,T}) where {T}
    # jacobian of vec(Q1'*Q2) wrt p2
    J = vcat(mrp_dcm_T_jacobian(p2,Q1[:,1]),
             mrp_dcm_T_jacobian(p2,Q1[:,2]),
             mrp_dcm_T_jacobian(p2,Q1[:,3]))
    idx = SA[1,4,7,2,5,8,3,6,9]
    return J[idx,:]
end

@inline function dcm_tilde_prod_jac_p1(p1::SVector{3,T},Q2::SMatrix{3,2,T,6}) where {T}
    # jacobian of vec(Q1'*Q2) wrt p1
    vcat(mrp_dcm_T_jacobian(p1,Q2[:,1])[SA[1,2],:],
         mrp_dcm_T_jacobian(p1,Q2[:,2])[SA[1,2],:])
end

@inline function dcm_tilde_prod_jac_p2(Q1::SMatrix{3,2,T,6},p2::SVector{3,T}) where {T}
    # jacobian of vec(Q1'*Q2) wrt p2
    J = vcat(mrp_dcm_T_jacobian(p2,Q1[:,1])[SA[1,2],:],
             mrp_dcm_T_jacobian(p2,Q1[:,2])[SA[1,2],:])
    idx = SA[1,3,2,4]
    return J[idx,:]
end
