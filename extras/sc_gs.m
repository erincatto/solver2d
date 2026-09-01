% This tests a two particle mass-spring-damper system with particle 1 connect to ground and
% particle 2 connected to particle 1. Particle 2 is much heavier than particle 1 and the
% simulation can go unstable if the soft constraint is made too stiff.
% Stability is improved by adding a relaxation step that applies the rigid constraint after
% the position update.
% Equilibrium at [0, -1]

% no relax stable up to 13.5 Hz
% with relax stable up to 20.5 Hz

relax = 0;

n = 2;

% Initially stretched: a vertical chain with rest spacing 1
ys = -(1:n)';
vs = zeros(n,1);

m = 1 * ones(n,1);
m(n) = 100;
invMasses = 1 ./ m;

% effective constraint masses: first constraint is particle1-ground, others between neighbours
km = zeros(n,1);
km(1) = invMasses(1);
km(2:n) = invMasses(1:n-1) + invMasses(2:n);
em = 1 ./ km;

h = 1/60;

% stiffness
hertz = 40;

% damping
zeta = 1;

% soft constraint parameters
omega = 2 * pi * hertz;
biasCoeff = omega / (2 * zeta + h * omega);
cc = h * omega * (2 * zeta + h * omega);
impulseCoeff = 1 / (1 + cc);
massCoeff = cc * impulseCoeff;

lambdas = zeros(n,1);
yys = ys;
cdot = zeros(n,1);
c = zeros(n,1);

for i = 1:1000
	% gravity/force
	vs += -10 * h;

	% warm start: apply existing constraint impulses to velocities
	% constraint 1: ground - particle 1
	vs(1:n-1) += invMasses(1:n-1) * (lambdas(1:n-1) - lambdas(2:n));
	vs(n) += invMasses(n) * lambdas(n);

	% sequential constraint solver
	c(1) = ys(1);
	c(2:n) = ys(2:n) - ys(1:n-1) + 1;

	for iter = 1 : 8
		cdot(1) = vs(1);
		dlambda = -massCoeff * em(1) * (cdot(1) + biasCoeff * c(1)) - impulseCoeff * lambdas(1);
		lambdas(1) += dlambda;
		vs(1) += invMasses(1) * dlambda;

		for j = 2:n
			cdot(j) = vs(j) - vs(j-1);
			dlambda = -massCoeff * em(j) * (cdot(j) + biasCoeff * c(j)) - impulseCoeff * lambdas(j);
			lambdas(j) += dlambda;

			vs(j - 1) -= invMasses(j - 1) * dlambda;
			vs(j) += invMasses(j) * dlambda;
		end
	end

	% integrate positions
	ys += h * vs;

	% relax (post-stabilization velocity correction)
	if relax == 1
		cdot(1) = vs(1);
		dlambda = -em(1) * cdot(1);
		lambdas(1) += dlambda;
		vs(1) += invMasses(1) * dlambda;

		for j = 2:n
			cdot(j) = vs(j) - vs(j-1);
			dlambda = -em(j) * cdot(j);
			lambdas(j) += dlambda;

			vs(j - 1) -= invMasses(j - 1) * dlambda;
			vs(j) += invMasses(j) * dlambda;
		end
	end

	yys = [yys ys];
end

plot(yys')
grid
