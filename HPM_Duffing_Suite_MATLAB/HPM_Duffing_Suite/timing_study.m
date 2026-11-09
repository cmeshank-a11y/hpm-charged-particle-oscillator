function Tm = timing_study(C, D, N, nrep)
%TIMING_STUDY  Wall-clock cost of each method on N points (median of nrep runs).
if nargin < 4, nrep = 50; end
t = linspace(0, 8, N); e = D/C; tau = sqrt(C)*t;
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
f = @(tt,y) [y(2); -C*y(1) - D*y(1)^3];
names = {}; times = [];
tt = zeros(1,nrep);
for r = 1:nrep, tic; ode45(f, t, [1;0], opts); tt(r) = toc; end
names{end+1} = 'ode45 (RelTol 1e-12)'; times(end+1) = median(tt);
for r = 1:nrep, tic; exact_cn(e, tau); tt(r) = toc; end
names{end+1} = 'Exact cn solution'; times(end+1) = median(tt);
for n = [0 1 2 4]
    for r = 1:nrep, tic; hpm_plain(e, tau, n); tt(r) = toc; end
    names{end+1} = sprintf('Plain HPM x^(%d)', n); times(end+1) = median(tt); %#ok<AGROW>
end
for n = [0 1 3]
    for r = 1:nrep, tic; hpm_fe(e, tau, n); tt(r) = toc; end
    names{end+1} = sprintf('HPM-FE order %d', n); times(end+1) = median(tt); %#ok<AGROW>
end
Tm.names = names; Tm.ms = 1e3*times; Tm.N = N; Tm.C = C; Tm.D = D;
end
