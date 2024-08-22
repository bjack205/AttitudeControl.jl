function euler_dynamics(x, u, t, p)
    q = Rotations.QuatRotation(@view x[1:4])
    h = x[5:7]
    τ = u
    J = p.J
    ω = J \ h
    hdot = τ - ω × h
    qdot = Rotations.kinematics(q, ω)    
    return [qdot; hdot]
end

attitude(x, p) = Rotations.QuatRotation(@view x[1:4])
momentum(x, p) = @view x[5:7]
angular_velocity(x, p) = p.J \ momentum(x)
inertial_momentum(x, p) = attitude(x, p) * momentum(x, p)