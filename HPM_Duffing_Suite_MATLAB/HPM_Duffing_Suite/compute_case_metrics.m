function R = compute_case_metrics(cases, T, N)
%COMPUTE_CASE_METRICS  Errors, energy drift and frequency errors for every case.
%  Reference = exact Jacobi-elliptic solution; an ode45 solution is also
%  computed (RelTol 1e-12, AbsTol 1e-14) as an independent numerical check.
if nargin < 2, T = 8; end
if nargin < 3, N = 501; end
t = linspace(0, T, N);
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
for i = 1:numel(cases)
    c = cases(i); e = c.eps; w = sqrt(c.C);
    tau = w*t;
    [Xe, dXe, Oex] = exact_cn(e, tau);
    Ee = energy_nd(e, Xe, dXe);
    % ODE reference in physical variables
    f = @(tt,y) [y(2); -c.C*y(1) - c.D*y(1)^3];
    [~, Y] = ode45(f, t, [c.A; 0], opts);
    Xo = Y(:,1).'/c.A; dXo = Y(:,2).'/(c.A*w);
    Eo = energy_nd(e, Xo, dXo);
    R(i).name = c.name; R(i).C = c.C; R(i).D = c.D; R(i).eps = e; %#ok<AGROW>
    R(i).periods = w*T/(2*pi);
    R(i).ode_maxdiff = max(abs(Xo - Xe));
    R(i).ode_dE = 100*(max(Eo)-min(Eo))/Eo(1);
    R(i).exact_dE = 100*(max(Ee)-min(Ee))/Ee(1);
    R(i).Omega_exact = Oex;
    % plain HPM, orders 0..4
    for n = 0:4
        [X, dX] = hpm_plain(e, tau, n);
        err = X - Xe;
        R(i).plain_MAE(n+1)  = max(abs(err));
        R(i).plain_RMSE(n+1) = sqrt(mean(err.^2));
        E = energy_nd(e, X, dX);
        R(i).plain_dE(n+1)   = 100*(max(E)-min(E))/E(1);
    end
    R(i).plain_ratio = R(i).plain_MAE(2:end)./R(i).plain_MAE(1:end-1);
    [R(i).plain_best_MAE, b] = min(R(i).plain_MAE); R(i).plain_best_order = b-1;
    % HPM with frequency expansion, orders 0..3
    for n = 0:3
        [X, dX, Om] = hpm_fe(e, tau, n);
        err = X - Xe;
        R(i).fe_MAE(n+1)  = max(abs(err));
        R(i).fe_RMSE(n+1) = sqrt(mean(err.^2));
        E = energy_nd(e, X, dX);
        R(i).fe_dE(n+1)   = 100*(max(E)-min(E))/E(1);
        R(i).fe_Omega_relerr(n+1) = 100*abs(Om/Oex - 1);
    end
    R(i).fe_ratio = R(i).fe_MAE(2:end)./R(i).fe_MAE(1:end-1);
    R(i).t = t; R(i).Xe = Xe; R(i).dXe = dXe; R(i).Xo = Xo; R(i).dXo = dXo;
end
end
