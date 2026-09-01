# This tests a two particle mass-spring-damper system with particle 1 connect to ground and
# particle 2 connected to particle 1. Particle 2 is much heavier than particle 1.
# Equilibrium at [0, -1]
# Direct solver

using LinearAlgebra
using Plots

# Initially stretched
ys = [-1.0; -2.0]
vs = [0.0; 0.0]
mA = 1.0
mB = 100.0
ms = [1/mA; 1/mB]
km = [ms[1]; ms[1] + ms[2]]
em = [1/km[1]; 1/km[2]]

# Jacobian
J = [1.0 0.0; -1.0 1.0]
invM = [ms[1] 0.0; 0.0 ms[2]]
K1 = J * invM * J'

h = 1/60

# stiffness
hertz = 30.0

# damping
zeta = 2.0

omega = 2 * pi * hertz
hw = h * omega
gammaA = km[1] / (hw * (2 * zeta + hw))
gammaB = km[2] / (hw * (2 * zeta + hw))
K2 = K1 + [gammaA 0.0; 0.0 gammaB]
Me = inv(K2)

gamma = [gammaA; gammaB]

beta = omega / (2 * zeta + hw)
cc = h * omega * (2 * zeta + h * omega)

lambdas = [0.0; 0.0]
yyd = copy(ys)

for i = 1:1000
    vs[1] += -10 * h
    vs[2] += -10 * h

    c = [ys[1]; ys[2] - ys[1] + 1]
    cdot = J * vs

    # Me includes the effect of gamma * lambda
    lambda = Me * (-cdot - beta * c)

    impulse = J' * lambda

    vs[1] += ms[1] * impulse[1]
    vs[2] += ms[2] * impulse[2]

    ys[1] += h * vs[1]
    ys[2] += h * vs[2]

    yyd = [yyd ys]
end

plot(yyd', title="direct")
