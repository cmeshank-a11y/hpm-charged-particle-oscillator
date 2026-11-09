function [X, dX, Omega, Tper] = exact_cn(eps_, tau)
%EXACT_CN  Exact solution of X'' + X + eps*X^3 = 0, X(0)=1, X'(0)=0 (eps >= 0).
%
%  X(tau) = cn(nu*tau | m),  nu = sqrt(1+eps),  m = eps/(2*(1+eps)),
%  Omega  = pi*nu/(2*K(m))   (angular frequency in units of sqrt(C)),
%  Tper   = 4*K(m)/nu        (period in tau).
%  See Byrd & Friedman (1971), Kovacic & Brennan (2011).
tau = tau(:).';
nu  = sqrt(1+eps_);
m   = eps_/(2*(1+eps_));
[sn, cn, dn] = ellipj(nu*tau, m);
X  = cn;
dX = -nu*sn.*dn;
K  = ellipke(m);
Omega = pi*nu/(2*K);
Tper  = 4*K/nu;
end
