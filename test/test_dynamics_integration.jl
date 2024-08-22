using AttitudeControl
using LinearAlgebra
using Rotations
using Rotations: params
using Test

J = Diagonal(rand(3))
q = rand(QuatRotation)
ω = rand(3)
h = J * ω
x = [params(q); h]
u = rand(3)
t = 0.0
p = (;J)
xdot = AttitudeControl.euler_dynamics(x, u, t, p)

# Zero momentum 
x = [params(q); h * 0]
xdot = AttitudeControl.euler_dynamics(x, u, t, p)
@test norm(xdot[1:4]) < 1e-8

# Zero torque
x = [params(q); h]
xdot = AttitudeControl.euler_dynamics(x, u * 0, t, p)
@test xdot[5:7] ≈ -ω × h

# Forward Euler
dt = 1/32
f(x,u,t) = AttitudeControl.euler_dynamics(x,u,t,p)
fd_rk4 = AttitudeControl.rk4(f, dt)

u = zeros(3)
fe_imid = AttitudeControl.implicit_midpoint(f, dt)
∇fe_imid = AttitudeControl.implicit_jacobians(fe_imid, dt) 
fd_imid = AttitudeControl.propagate_implicit(fe_imid, ∇fe_imid[2], tol=1e-16)
x2 = fd_imid(x, u, t)
x2 = fd_rk4(x, u, t)

hI_1 = AttitudeControl.inertial_momentum(x, p)
hI_2 = AttitudeControl.inertial_momentum(x2, p)
hI_2 - hI_1

x0 = [params(q); h]
u = zeros(3)
x = copy(x0)
for i = 1:1000
    x = fd_imid(x, u, t)
    # x = fd_rk4(x, u, t)
end
hI_0 = AttitudeControl.inertial_momentum(x0, p)
hI_2 = AttitudeControl.inertial_momentum(x, p)
hI_2 - hI_0
norm(x[1:4])