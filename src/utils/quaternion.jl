function dcm_from_q(q::SVector{4,T}) where {T}
    #DCM from quaternion, hamilton product, scalar first
    # pull our the parameters from the quaternion
    q4,q1,q2,q3 = normalize(q)

    # DCM
    Q = @SArray [(2*q1^2+2*q4^2-1)   2*(q1*q2 - q3*q4)   2*(q1*q3 + q2*q4);
          2*(q1*q2 + q3*q4)  (2*q2^2+2*q4^2-1)   2*(q2*q3 - q1*q4);
          2*(q1*q3 - q2*q4)   2*(q2*q3 + q1*q4)  (2*q3^2+2*q4^2-1)]
end

function dcm_tilde_from_q(q::SVector{4,T}) where {T}
    #DCM from quaternion, hamilton product, scalar first
    # pull our the parameters from the quaternion
    q4,q1,q2,q3 = normalize(q)

    # DCM
    Q = @SArray [(2*q1^2+2*q4^2-1)   2*(q1*q2 - q3*q4) ;
          2*(q1*q2 + q3*q4)  (2*q2^2+2*q4^2-1)   ;
          2*(q1*q3 - q2*q4)   2*(q2*q3 + q1*q4)  ]
end


function p_from_q(q::SVector{4,T}) where {T}
    return q[SA[2,3,4]]/(1+q[1])
end

function q_from_p(p::SVector{3,T}) where {T}
    return (1/(1+dot(p,p)))*vcat((1-dot(p,p)),2*p)
end
