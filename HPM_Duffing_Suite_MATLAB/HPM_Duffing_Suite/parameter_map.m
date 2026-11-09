function M = parameter_map(neps, nper, Npts)
%PARAMETER_MAP  Error maps in the scaled plane (eps, number of linear periods).
%  Because tau = sqrt(C)*t and eps = D*A^2/C, the accuracy of every method
%  over a window t in [0,T] depends only on eps and n_p = sqrt(C)*T/(2*pi).
if nargin < 1, neps = 81; end
if nargin < 2, nper = 60; end
if nargin < 3, Npts = 1201; end
M.eps = logspace(-3, 1, neps);
M.np  = linspace(0.25, 10, nper);
M.plain_MAE = zeros(nper, neps, 5);
M.fe_MAE    = zeros(nper, neps, 4);
for j = 1:neps
    e = M.eps(j);
    for i = 1:nper
        tau = linspace(0, 2*pi*M.np(i), Npts);
        Xe = exact_cn(e, tau);
        [~,~,Tm] = hpm_plain(e, tau, 4);
        S = cumsum(Tm, 1);
        for n = 0:4, M.plain_MAE(i,j,n+1) = max(abs(S(n+1,:) - Xe)); end
        for n = 0:3, M.fe_MAE(i,j,n+1) = max(abs(hpm_fe(e, tau, n) - Xe)); end
    end
end
[~, b] = min(M.plain_MAE, [], 3); M.plain_best = b - 1;
[~, b] = min(M.fe_MAE,    [], 3); M.fe_best    = b - 1;
end
