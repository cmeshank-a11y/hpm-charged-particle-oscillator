function write_tables(R, Tm, outdir)
%WRITE_TABLES  Writes the manuscript tables as CSV files.
w = @(name) fopen(fullfile(outdir, name), 'w');
f = w('Table1_plain_HPM_MAE.csv');
fprintf(f, 'case,C,D,eps,periods,MAE_x0,MAE_x1,MAE_x2,MAE_x3,MAE_x4,improvement_x0_over_x4,best_order\n');
for r = R, fprintf(f, '"%s",%g,%g,%g,%.4f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%d\n', r.name, r.C, r.D, r.eps, r.periods, r.plain_MAE, r.plain_MAE(1)/r.plain_MAE(5), r.plain_best_order); end
fclose(f);
f = w('Table2_plain_HPM_RMSE.csv');
fprintf(f, 'case,C,D,RMSE_x0,RMSE_x1,RMSE_x2,RMSE_x3,RMSE_x4\n');
for r = R, fprintf(f, '"%s",%g,%g,%.6e,%.6e,%.6e,%.6e,%.6e\n', r.name, r.C, r.D, r.plain_RMSE); end
fclose(f);
f = w('Table3_HPM_FE_MAE.csv');
fprintf(f, 'case,C,D,MAE_FE0,MAE_FE1,MAE_FE2,MAE_FE3,Omega_exact,Omega_err_pct_FE0,Omega_err_pct_FE1,Omega_err_pct_FE2,Omega_err_pct_FE3\n');
for r = R, fprintf(f, '"%s",%g,%g,%.6e,%.6e,%.6e,%.6e,%.8f,%.6e,%.6e,%.6e,%.6e\n', r.name, r.C, r.D, r.fe_MAE, r.Omega_exact, r.fe_Omega_relerr); end
fclose(f);
f = w('Table4_energy.csv');
fprintf(f, 'case,C,D,dE_exact_pct,dE_ode45_pct,dE_plain_x1_pct,dE_plain_x4_pct,dE_FE0_pct,dE_FE3_pct,ode45_vs_exact_maxdiff\n');
for r = R, fprintf(f, '"%s",%g,%g,%.3e,%.3e,%.6e,%.6e,%.6e,%.6e,%.3e\n', r.name, r.C, r.D, r.exact_dE, r.ode_dE, r.plain_dE(2), r.plain_dE(5), r.fe_dE(1), r.fe_dE(4), r.ode_maxdiff); end
fclose(f);
f = w('Table5_error_ratios.csv');
fprintf(f, 'case,C,D,plain_r1,plain_r2,plain_r3,plain_r4,fe_r1,fe_r2,fe_r3\n');
for r = R, fprintf(f, '"%s",%g,%g,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', r.name, r.C, r.D, r.plain_ratio, r.fe_ratio); end
fclose(f);
f = w('Table6_timing.csv');
fprintf(f, 'method,time_ms\n');
for k = 1:numel(Tm.names), fprintf(f, '%s,%.6f\n', Tm.names{k}, Tm.ms(k)); end
fclose(f);
end
