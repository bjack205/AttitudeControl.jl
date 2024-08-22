function simulate(discrete_dynamics::Function, ctrl::Function, x0; tspan=(0,1), dt=0.01)
    u0 = ctrl(x0, tspan[1])
    n = length(x0)
    m = length(u0)

    t = range(tspan[1], tspan[2], step=dt)
    N = length(t) 
    x_traj = [zeros(n) for k = 1:N]
    u_traj = [zeros(m) for k = 1:N-1]
    x_traj[1] .= x0
    for k = 1:N-1
        u[k] = ctrl(x[k], t[k])
        x[k+1] = discrete_dynamics(x[k], u[k], t[k])
    end
    return t,x,u
end
