% This tests a two particle mass-spring-damper system with particle 1 connect to ground and 
% particle 2 connected to particle 1. Particle 2 is much heavier than particle 1 and the
% simulation can go unstable if the soft constraint is made too stiff.
% Stability is improved by adding a relaxation step that applies the rigid constraint after
% the position update.

% no relax stable up to 13.5 Hz
% with relax stable up to 20.5 Hz
relax = 1;
hertz = 20.5;

ys = [0;-1];
vs = [0;0];
ms = [0.167;0.000995];
km = [ms(1);ms(1) + ms(2)];
em = [1/km(1); 1/km(2)];
h = 1/60;
% no damping
zeta = 0;
omega = 2 * pi * hertz;
biasCoeff = omega / (2 * zeta + h * omega);
c = h * omega * (2 * zeta + h * omega);
impulseCoeff = 1 / (1 + c);
massCoeff = c * impulseCoeff;

lambdas = [0;0];
yys = [ys];

for i = 1:10000
	vs(1) += -10 * h;
	vs(2) += -10 * h;

	% warm start
	vs(1) += ms(1) * (lambdas(1) - lambdas(2));
	vs(2) += ms(2) * lambdas(2);

	cdot1 = vs(1);
	c1 = ys(1);
	bias = biasCoeff * c1;
	lambda1 = -massCoeff * em(1) * (cdot1 + bias) - impulseCoeff * lambdas(1);
	lambdas(1) += lambda1;

	vs(1) += ms(1) * lambda1;

	cdot2 = vs(2) - vs(1);
	c2 = ys(2) - ys(1) + 1;
	lambda2 = -massCoeff * em(2) * (cdot2 + biasCoeff * c2) - impulseCoeff * lambdas(2);
	lambdas(2) += lambda2;

	vs(1) -= ms(1) * lambda2;
	vs(2) += ms(2) * lambda2;

	ys(1) += h * vs(1);
	ys(2) += h * vs(2);

	% relax
	if relax == 1
		cdot1 = vs(1);
		lambda1 = -em(1) * cdot1;
		lambdas(1) += lambda1;

		vs(1) += ms(1) * lambda1;

		cdot2 = vs(2) - vs(1);
		lambda2 = -em(2) * cdot2;
		lambdas(2) += lambda2;

		vs(1) -= ms(1) * lambda2;
		vs(2) += ms(2) * lambda2;
	end

	yys = [yys ys];
end

plot(yys')

