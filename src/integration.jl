function rk4(f,dt)
    return function (x,u,t)
        k1 = f(x,u,t) * dt
        k2 = f(x + k1 / 2, u, t + dt / 2) * dt
        k3 = f(x + k2 / 2, u, t + dt / 2) * dt
        k4 = f(x + k3, u, t + dt) * dt
        x + (k1 + 2k2 + 2k3 + k4) / 6
    end
end

function forward_euler(f, dt)
    return function (x,u,t)
        x + f(x,u,t) * dt
    end
end

function explicit_midpoint(f, dt)
    return function (x,u,t)
        x + f(x + f(x,u,t) * dt / 2, u, t + dt / 2) * dt
    end
end

function implicit_midpoint(f, dt)
    return function (x1,x2,u,t)
        x_mid = (x1 + x2) / 2
        x1 + f(x_mid,u,t + dt / 2) * dt - x2
    end
end

function explicit_jacobians(f, dt)
    dx(x,u,t) = ForwardDiff.jacobian(x->f(x,u,t), x)
    du(x,u,t) = ForwardDiff.jacobian(u->f(x,u,t), u)
    return dx, du
end

function implicit_jacobians(fe, dt)
    dx1(x1,x2,u,t) = ForwardDiff.jacobian(x1->fe(x1,x2,u,t), x1)
    dx2(x1,x2,u,t) = ForwardDiff.jacobian(x2->fe(x1,x2,u,t), x2)
    du(x1,x2,u,t) = ForwardDiff.jacobian(u->fe(x1,x2,u,t), u)
    return dx1, dx2, du
end

function propagate_implicit(fd, dfx2; max_iters=20, tol=1e-6)
    return function (x,u,t)
        x1 = copy(x)
        x2 = copy(x) 
        did_succeed = false
        for iter = 1:max_iters
            r = fd(x1,x2,u,t)
            @show norm(r)
            if norm(r) < tol
                did_succeed = true
                return x2
            end
            ∇r = dfx2(x1, x2, u, t)
            dx = -∇r \ r
            x2 .+= dx
        end
        if !did_succeed
            @warn("propagate_implicit did not converge")
        end
        return x2
    end
end