function make_figures(R, M, figdir)
%MAKE_FIGURES  Produces Figures 1-5 of the manuscript as 300-dpi PNG files.
%  Fig1 solutions/errors/phase, Fig2 convergence, Fig3 energy, Fig4 maps, Fig5 frequency.
if nargin < 3, figdir = 'figures'; end
if ~exist(figdir, 'dir'), mkdir(figdir); end
set(0, 'defaultaxesfontsize', 9, 'defaulttextfontsize', 9, 'defaultlinelinewidth', 1.1);
col = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188; 0.494 0.184 0.556; 0.929 0.694 0.125];
fig1(R, col, figdir); fig2(R, col, figdir); fig3(R, col, figdir);
fig4(R, M, figdir);   fig5(col, figdir);
end
% -------------------------------------------------------------------------
function savefig_(f, name, w, h, figdir)
set(f, 'PaperUnits', 'inches', 'PaperPosition', [0 0 w h], 'PaperSize', [w h]);
print(f, fullfile(figdir, [name '.png']), '-dpng', '-r300');
close(f);
end
% -------------------------------------------------------------------------
function logticks_(ax, axname, ex)
% clean 10^k tick labels (avoids zero-padded exponents in some renderers)
lab = arrayfun(@(k) sprintf('10^{%d}', k), ex, 'uniformoutput', false);
set(ax, [axname 'tick'], 10.^ex, [axname 'ticklabel'], lab);
end
% -------------------------------------------------------------------------
function fig1(R, col, figdir)
% time histories, error growth and phase portraits for three regimes
idx = [1 3 7]; lab = 'abcdefghi';
f = figure('visible', 'off');
for r = 1:3
    c = R(idx(r)); t = c.t; tau = sqrt(c.C)*t; e = c.eps;
    [X1]      = hpm_plain(e, tau, 1);
    [X4, dX4] = hpm_plain(e, tau, 4);
    [F3, dF3] = hpm_fe(e, tau, 3);
    [F0]      = hpm_fe(e, tau, 0);
    % (1) solutions
    subplot(3,3,3*r-2);
    plot(t, c.Xe, 'k-', 'linewidth', 1.6); hold on;
    plot(t, X1, '--', 'color', col(1,:)); plot(t, X4, '-', 'color', col(2,:));
    plot(t, F3, '-.', 'color', col(3,:));
    ylim([-1.6 1.6]); xlim([0 8]); grid on; xlabel('t'); ylabel('x(t)');
    title(sprintf('(%s) %s, \\epsilon = %.3g', lab(3*r-2), c.name, e));
    % (2) absolute errors
    subplot(3,3,3*r-1);
    [~,~,Tm] = hpm_plain(e, tau, 4); S = cumsum(Tm,1);
    for n = 0:4
        semilogy(t, abs(S(n+1,:) - c.Xe) + 1e-16, '-', 'color', col(n+1,:)); hold on;
    end
    semilogy(t, abs(F0 - c.Xe) + 1e-16, 'k:', 'linewidth', 1.3);
    semilogy(t, abs(F3 - c.Xe) + 1e-16, 'k-', 'linewidth', 1.3);
    xlim([0 8]); ylim([1e-8 1e3]); grid on; xlabel('t'); ylabel('|x - x_{exact}|'); logticks_(gca, 'y', -8:2:2);
    title(sprintf('(%s) Absolute error', lab(3*r-1)));
    % (3) phase portraits (non-dimensional velocity dx/dtau)
    subplot(3,3,3*r);
    plot(c.Xe, c.dXe, 'k-', 'linewidth', 1.6); hold on;
    plot(X4, dX4, '-', 'color', col(2,:)); plot(F3, dF3, '-.', 'color', col(3,:));
    axis([-2 2 -3 3]); grid on; xlabel('x'); ylabel('dx/d\tau');
    title(sprintf('(%s) Phase plane', lab(3*r)));
end
savefig_(f, 'Fig1_solutions_errors_phase', 7.2, 7.6, figdir);
end
% -------------------------------------------------------------------------
function fig2(R, col, figdir)
% relative energy deviation |E(t)-E(0)|/E(0)
idx = [1 3 5 7]; lab = 'abcd';
f = figure('visible', 'off');
for k = 1:4
    c = R(idx(k)); t = c.t; tau = sqrt(c.C)*t; e = c.eps;
    Eo = energy_nd(e, c.Xo, c.dXo);
    [X4, dX4] = hpm_plain(e, tau, 4); E4 = energy_nd(e, X4, dX4);
    [F3, dF3] = hpm_fe(e, tau, 3);    E3 = energy_nd(e, F3, dF3);
    [F0, dF0] = hpm_fe(e, tau, 0);    E0 = energy_nd(e, F0, dF0);
    subplot(2,2,k);
    semilogy(t, abs(Eo/Eo(1) - 1) + 1e-17, 'k-'); hold on;
    semilogy(t, abs(E4/E4(1) - 1) + 1e-17, '-', 'color', col(2,:));
    semilogy(t, abs(E0/E0(1) - 1) + 1e-17, ':', 'color', col(1,:), 'linewidth', 1.4);
    semilogy(t, abs(E3/E3(1) - 1) + 1e-17, '-.', 'color', col(3,:));
    xlim([0 8]); ylim([1e-16 1e10]); grid on; xlabel('t'); ylabel('|E(t)/E(0) - 1|');
    title(sprintf('(%s) %s, \\epsilon = %.3g', lab(k), c.name, e));
    logticks_(gca, 'y', -15:5:10);
    if k == 1, legend({'ode45','HPM x^{(4)}','HPM-FE (0)','HPM-FE (3)'}, 'location', 'north', 'orientation', 'horizontal', 'numcolumns', 2, 'fontsize', 7); end
end
savefig_(f, 'Fig3_energy', 7.0, 5.4, figdir);
end
% -------------------------------------------------------------------------
function fig3(R, col, figdir)
% error versus truncation order for all cases
mk = 'osd^v<>';
cc = [col; 0.301 0.745 0.933; 0.635 0.078 0.184];
f = figure('visible', 'off');
ax1 = subplot(1,2,1);
for i = 1:numel(R)
    semilogy(0:4, R(i).plain_MAE, ['-' mk(i)], 'color', cc(i,:), 'markerfacecolor', cc(i,:), 'markersize', 4); hold on;
end
grid on; xlabel('Truncation order n'); ylabel('Maximum error'); xlim([-0.2 4.2]); set(gca, 'xtick', 0:4);
title('(a) Plain HPM'); ylim([1e-5 1e3]); logticks_(gca, 'y', -5:3);
ax2 = subplot(1,2,2);
for i = 1:numel(R)
    semilogy(0:3, R(i).fe_MAE, ['-' mk(i)], 'color', cc(i,:), 'markerfacecolor', cc(i,:), 'markersize', 4); hold on;
end
grid on; xlabel('Truncation order n'); ylabel('Maximum error'); xlim([-0.2 3.2]); set(gca, 'xtick', 0:3);
title('(b) HPM with frequency expansion'); ylim([1e-7 1e1]); logticks_(gca, 'y', -7:1);
L = arrayfun(@(r) sprintf('%s (\\epsilon=%.3g)', r.name, r.eps), R, 'uniformoutput', false);
lg = legend(L, 'location', 'eastoutside', 'fontsize', 7);
set(ax1, 'position', [0.07 0.15 0.31 0.75]); set(ax2, 'position', [0.46 0.15 0.31 0.75]);
savefig_(f, 'Fig2_convergence', 9.5, 3.6, figdir);
end
% -------------------------------------------------------------------------
function fig4(R, M, figdir)
% maps in the (log10 eps, n_p) plane
le = log10(M.eps); np = M.np;
f = figure('visible', 'off');
A = log10(M.plain_MAE(:,:,5));
B = M.plain_best;
G = log10(M.plain_MAE(:,:,5) ./ M.plain_MAE(:,:,1));
Cf = log10(M.fe_MAE(:,:,4));
pts = [log10([R.eps]); [R.periods]];
subplot(1,3,1);
imagesc(le, np, A); axis xy; caxis([-8 3]); colorbar; hold on;
contour(le, np, G, [0 0], 'w-', 'linewidth', 1.6);
contour(le, np, A, [-2 -2], 'k--', 'linewidth', 1.0);
plot(pts(1,:), pts(2,:), 'wo', 'markerfacecolor', 'k', 'markersize', 4);
xlabel('log_{10} \epsilon'); ylabel('Number of periods n_p'); title('(a) Plain HPM x^{(4)}: log_{10} MAE');
subplot(1,3,2);
imagesc(le, np, B); axis xy; caxis([-0.5 4.5]); colormap(gca, [0.85 0.33 0.10; 0.93 0.69 0.13; 0.47 0.67 0.19; 0.30 0.75 0.93; 0.00 0.45 0.74]); h = colorbar; set(h, 'ytick', 0:4);
hold on; plot(pts(1,:), pts(2,:), 'wo', 'markerfacecolor', 'k', 'markersize', 4);
xlabel('log_{10} \epsilon'); ylabel('Number of periods n_p'); title('(b) Plain HPM: optimal order');
subplot(1,3,3);
imagesc(le, np, Cf); axis xy; caxis([-8 3]); colorbar; hold on;
contour(le, np, Cf, [-2 -2], 'k--', 'linewidth', 1.0);
plot(pts(1,:), pts(2,:), 'wo', 'markerfacecolor', 'k', 'markersize', 4);
xlabel('log_{10} \epsilon'); ylabel('Number of periods n_p'); title('(c) HPM-FE order 3: log_{10} MAE');
savefig_(f, 'Fig4_parameter_maps', 10.5, 3.9, figdir);
end
% -------------------------------------------------------------------------
function fig5(col, figdir)
% frequency accuracy and long-window amplitude accuracy versus eps
e = logspace(-3, 1, 200);
Fe = nan(4, numel(e)); Ae = nan(5, numel(e));
tau = linspace(0, 2*pi*5, 2001);       % five linear periods
for j = 1:numel(e)
    [Xe,~,Oex] = exact_cn(e(j), tau);
    for n = 0:3
        [X,~,Om] = hpm_fe(e(j), tau, n);
        Fe(n+1,j) = abs(Om/Oex - 1);
        Ae(n+1,j) = max(abs(X - Xe));
    end
    Ae(5,j) = max(abs(hpm_plain(e(j), tau, 4) - Xe));
end
f = figure('visible', 'off');
subplot(1,2,1);
for n = 0:3, loglog(e, Fe(n+1,:), '-', 'color', col(n+1,:)); hold on; end
plot([1 1], [1e-16 10], 'k:');
grid on; xlabel('\epsilon'); ylabel('|\Omega_{FE}/\Omega_{exact} - 1|'); ylim([1e-14 10]); logticks_(gca, 'y', -14:2:0); logticks_(gca, 'x', -3:1);
legend({'order 0','order 1','order 2','order 3'}, 'location', 'southeast', 'fontsize', 7);
title('(a) Frequency error, HPM-FE');
subplot(1,2,2);
for n = 0:3, loglog(e, Ae(n+1,:), '-', 'color', col(n+1,:)); hold on; end
loglog(e, Ae(5,:), 'k--');
plot([1 1], [1e-16 1e4], 'k:');
grid on; xlabel('\epsilon'); ylabel('Maximum error over 5 periods'); ylim([1e-12 1e4]); logticks_(gca, 'y', -12:2:4); logticks_(gca, 'x', -3:1);
legend({'FE order 0','FE order 1','FE order 2','FE order 3','plain HPM x^{(4)}'}, 'location', 'southeast', 'fontsize', 7);
title('(b) Amplitude error over n_p = 5');
savefig_(f, 'Fig5_frequency', 7.2, 3.3, figdir);
end
