function [V, dV] = hpm_hierarchy_ode(tau)
%HPM_HIERARCHY_ODE  Numerically integrates the HPM recurrences (eps factored out)
%   u0''+u0 = 0,           u0(0)=1
%   u1''+u1 = -u0^3
%   u2''+u2 = -3u0^2u1
%   u3''+u3 = -3(u0^2u2 + u0u1^2)
%   u4''+u4 = -(3u0^2u3 + 6u0u1u2 + u1^3),   un(0)=un'(0)=0 for n>=1.
%  Used only to verify the closed-form coefficients in hpm_coefficients.m.
f = @(s,y) rhs(y);
y0 = zeros(10,1); y0(1) = 1;
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
[~, Y] = ode45(f, tau(:), y0, opts);
V  = Y(:,1:2:end).';  dV = Y(:,2:2:end).';
end
function dy = rhs(y)
u = y(1:2:end); du = y(2:2:end);
F = [0; -u(1)^3; -3*u(1)^2*u(2); -3*(u(1)^2*u(3)+u(1)*u(2)^2); ...
     -(3*u(1)^2*u(4) + 6*u(1)*u(2)*u(3) + u(2)^3)];
dy = zeros(10,1);
dy(1:2:end) = du;
dy(2:2:end) = -u + F;
end
