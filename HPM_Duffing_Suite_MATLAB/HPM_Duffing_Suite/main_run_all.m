%MAIN_RUN_ALL  Reproduces every number, table and figure of the manuscript
%   "Accuracy, convergence and failure of the homotopy perturbation method
%    for the nonlinear charged particle oscillation equation".
%
%   Runs in MATLAB (R2019b or later) and GNU Octave (>= 7).
%   Outputs:  results/*.csv, results/summary.txt, results/results.mat,
%             figures/Fig1..Fig5 (*.png, 300 dpi).
clear; clc;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
cd(here);
if ~exist('results','dir'), mkdir('results'); end
if ~exist('figures','dir'), mkdir('figures'); end
fid = fopen(fullfile('results','summary.txt'), 'w');
fprintf(fid, 'HPM for x'''' + C x + D x^3 = 0, x(0)=1, x''(0)=0.  Run on %s\n\n', datestr(now));

%% 1. Verification of the analytical expressions
V = verify_terms(fid); verify_terms(1);

%% 2. Seven parameter sets, window t in [0,8], N = 501
cases = case_definitions();
R = compute_case_metrics(cases, 8, 501);

%% 3. Scaling check: same eps and same number of periods -> identical errors
c1 = R(5);                                   % C = 1, D = 1, T = 8
cs = struct('name','scaled','C',0.1,'D',0.1,'eps',1,'A',1);
Rs = compute_case_metrics(cs, 8*sqrt(1/0.1), 501);
fprintf(fid, '\nScaling check (C,D,T)=(1,1,8) vs (0.1,0.1,%.4f): max diff in plain MAE = %.2e, in FE MAE = %.2e\n', ...
    8*sqrt(10), max(abs(c1.plain_MAE - Rs.plain_MAE)), max(abs(c1.fe_MAE - Rs.fe_MAE)));

%% 4. Parameter map in the scaled plane (eps, n_p)
M = parameter_map(81, 60, 1201);
G = M.plain_MAE(:,:,5) < M.plain_MAE(:,:,1);
fprintf(fid, 'Parameter map: plain x^(4) beats x^(0) on %.1f%% of the (eps,n_p) grid\n', 100*mean(G(:)));
for np_ = [1 2 5 10]
    [~, i] = min(abs(M.np - np_));
    jb = find(M.plain_MAE(i,:,5) >= M.plain_MAE(i,:,1), 1, 'first');
    j2 = find(M.plain_MAE(i,:,5) > 1e-2, 1, 'first');
    j3 = find(M.fe_MAE(i,:,4) > 1e-2, 1, 'first');
    fprintf(fid, '  n_p = %5.2f : plain x4 stops improving on x0 at eps = %.3g ; plain x4 MAE>1e-2 at eps = %.3g ; FE3 MAE>1e-2 at eps = %.3g\n', ...
        M.np(i), M.eps(jb), M.eps(j2), M.eps(j3));
end
% best and worst plain HPM points of the map
A = M.plain_MAE(:,:,5);
[mn, k] = min(A(:)); [i,j] = ind2sub(size(A), k);
fprintf(fid, '  best plain x4: MAE=%.2e at eps=%.3g, n_p=%.2f\n', mn, M.eps(j), M.np(i));
[mx, k] = max(A(:)); [i,j] = ind2sub(size(A), k);
fprintf(fid, '  worst plain x4: MAE=%.2e at eps=%.3g, n_p=%.2f\n', mx, M.eps(j), M.np(i));
F = M.fe_MAE(:,:,4); [mx, k] = max(F(:)); [i,j] = ind2sub(size(F), k);
fprintf(fid, '  worst FE3: MAE=%.2e at eps=%.3g, n_p=%.2f\n', mx, M.eps(j), M.np(i));
F0 = M.fe_MAE(:,:,1); fprintf(fid, '  worst FE0 over whole map: MAE=%.3f\n', max(F0(:)));

%% 5. Frequency of HPM-FE for large eps
fprintf(fid, '\nHPM-FE frequency relative error (%%) for eps = 1, 2, 4, 10, 100:\n');
for e = [1 2 4 10 100]
    [~,~,Oex] = exact_cn(e, 0); r = nan(1,4);
    for n = 0:3, [~,~,Om] = hpm_fe(e, 0, n); r(n+1) = 100*abs(Om/Oex-1); end
    fprintf(fid, '  eps=%6g : %s\n', e, mat2str(r, 4));
end
[~,~,Oinf] = exact_cn(1e8, 0);
fprintf(fid, '  eps->inf limit of order-0 (He) frequency error: %.3f %%\n', 100*abs(sqrt(3/4)*sqrt(1e8)/Oinf - 1));

%% 6. Computational cost
Tm = timing_study(2.0, 1.0, 1000, 50);
fprintf(fid, '\nTiming (N = %d, C = %.1f, D = %.1f), median of 50 runs:\n', Tm.N, Tm.C, Tm.D);
for k = 1:numel(Tm.names), fprintf(fid, '  %-24s %9.3f ms\n', Tm.names{k}, Tm.ms(k)); end

%% 7. Write tables
write_tables(R, Tm, 'results');
fprintf(fid, '\nPer-case summary:\n');
for i = 1:numel(R)
    fprintf(fid, '%-16s eps=%.3f n_p=%.2f | ode45-exact %.1e | best plain order %d (MAE %.2e) | best FE order MAE %.2e\n', ...
        R(i).name, R(i).eps, R(i).periods, R(i).ode_maxdiff, R(i).plain_best_order, R(i).plain_best_MAE, min(R(i).fe_MAE));
end
fclose(fid);
Rlite = rmfield(R, {'t','Xe','dXe','Xo','dXo'});
save(fullfile('results','results.mat'), 'Rlite', 'M', 'Tm', 'V');

%% 8. Figures
make_figures(R, M, 'figures');
type(fullfile('results','summary.txt'));
