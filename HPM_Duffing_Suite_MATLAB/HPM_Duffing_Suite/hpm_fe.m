function [X, dX, Omega] = hpm_fe(eps_, tau, order)
%HPM_FE  HPM with frequency expansion (He's parameter-expanded HPM).
%
%  The frequency is embedded in the homotopy, Omega^2 = 1 + p*w1 + p^2*w2 + ...,
%  and each w_j is fixed by removing the secular (resonant) term at order p^j.
%  The order-n approximation is
%      X = sum_{k=0}^{n} eps^k X_k(Omega_n*tau),
%      Omega_n^2 = 1 + sum_{j=1}^{n+1} w_j,
%  i.e. the frequency is always one order ahead of the amplitude.  order <= 3.
%  Omega is returned in units of the linear frequency sqrt(C).
if order > 3, error('hpm_fe: order must be <= 3'); end
[~, XF, wF] = hpm_coefficients();
tau = tau(:).';
Om2 = 1 + sum(wF(1:order+1) .* eps_.^(1:order+1));
if Om2 <= 0
    X = nan(size(tau)); dX = X; Omega = NaN; return;   % truncated series has failed
end
Omega = sqrt(Om2);
s  = Omega*tau;
X  = zeros(size(tau)); dXs = zeros(size(tau));
for r = 1:size(XF,1)
    n = XF(r,1); if n > order, continue; end
    c = XF(r,2)/XF(r,3); m = XF(r,4);
    X   = X   + eps_^n * c * cos(m*s);
    dXs = dXs - eps_^n * c * m * sin(m*s);
end
dX = Omega*dXs;
end
