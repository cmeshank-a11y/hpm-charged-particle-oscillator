function out = verify_terms(fid)
%VERIFY_TERMS  Independent checks of the analytical HPM expressions.
%  1) closed-form u_n vs. direct numerical integration of the recurrences;
%  2) initial conditions u_n(0)=u_n'(0)=0 (n>=1);
%  3) HPM-FE frequency series vs. exact elliptic frequency for small eps.
if nargin < 1, fid = 1; end
tau = linspace(0, 8, 801);
[V, dV] = hpm_hierarchy_ode(tau);
U = hpm_coefficients(); %#ok<NASGU>
out.maxdiff = zeros(1,5);
for n = 0:4
    [~,~,T]  = hpm_plain(1, tau, 4);     % eps = 1 isolates u_n
    out.maxdiff(n+1) = max(abs(T(n+1,:) - V(n+1,:)));
end
fprintf(fid, 'Check 1: max |u_n(closed form) - u_n(ODE)| on tau in [0,8]:\n');
fprintf(fid, '   n=%d : %.2e\n', [0:4; out.maxdiff]);
[~,~,T0] = hpm_plain(1, 0, 4);
dt = 1e-5; [~,~,Tp] = hpm_plain(1, dt, 4); [~,~,Tn] = hpm_plain(1, -dt, 4);
fprintf(fid, 'Check 2: u_n(0) = %s ;  central-difference u_n''(0) = %s\n', ...
    mat2str(T0(:).',3), mat2str(((Tp(:)-Tn(:))/(2*dt)).',3));
e = [1e-3 1e-2 1e-1];
fprintf(fid, 'Check 3: HPM-FE frequency vs exact (relative error):\n');
for k = 1:numel(e)
    [~,~,Oex] = exact_cn(e(k), 0);
    r = zeros(1,4);
    for n = 0:3, [~,~,Om] = hpm_fe(e(k), 0, n); r(n+1) = abs(Om/Oex-1); end
    fprintf(fid, '   eps=%g : orders 0..3 -> %s\n', e(k), mat2str(r,3));
end
end
