clear
close all
clc
% parameter
R2 = 6371e3;
R1 = 6361e3;
M = 5.965e24;
a = R2;
x = R1 : R2;
RE = @(x, xref) log10(abs(x ./ xref - 1));
sqrtSigma = @(x) log10( sqrt(degVar(x)) );
r0 = R1;
cmp = [0.85, 0.10, 0.10; 0.00, 0.45, 0.74; 0.20, 0.60, 0.20];
%% poly
rho_poly = [3000, -1.7, 4.1 * 1e-4, -4.1 * 1e-8, 1.42 * 1e-12];
frho_poly = @(xr) rho_poly(1) + rho_poly(2) * (xr - r0) + rho_poly(3) * (xr - r0) .^ 2 + rho_poly(4) * (xr - r0) .^ 3 + rho_poly(5) * (xr - r0) .^ 4;
vnm_poly_ref = shell_C00_poly(R1, R2, rho_poly, M, r0);
%% exp
rho0 = 3 * 1e3;
mu = 7 * 1e-4;
frho_exp = @(xr) rho0 * exp(-mu * (xr - r0));
rho_exp = [rho0, mu];
vnm_exp_ref = shell_C00_exp(R1, R2, rho_exp, M, r0);
%% para
rho0 = 3000;
mu = -1.4;
frho_par = @(xr) rho0 ^ 3 ./ (rho0 - mu * (xr - r0)) .^ 2;
rho_par = [rho0, mu];
vnm_par_ref = shell_C00_para(R1, R2, rho_par, M, r0);
%% plot
fontsize = 9;
linewidth = 1.2;
fig = figure;
tiledlayout(1, 1, "Padding", "tight", 'TileSpacing', 'tight', Parent=fig);
nexttile
plot(x * 1e-3, frho_poly(x'), "LineWidth", linewidth);
hold on
plot(x * 1e-3, frho_exp(x), "LineWidth", linewidth);
plot(x * 1e-3, frho_par(x), "LineWidth", linewidth);
xlim([R1, R2] * 1e-3);
xticks(6361 : 2 : 6371)
xlabel('Radius (km)', 'Interpreter', 'latex');
ylabel('Radial density functions ($\mathrm{kg~m^{-3}}$)', 'Interpreter', 'latex');
grid on
legend('Polynomial', 'Exponential', 'Parabolic', 'fontsize', fontsize, 'NumColumns', 1, 'Interpreter', 'latex');
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);
set(gcf,'Units','centimeters','Position', [1.3, 3, 8, 6]);
exportgraphics(gcf, '../manuscript/marked/density.pdf');
% dis tesseroid
dlat = 1;
nmax = 180;
lon0 = -180 + dlat / 2 : dlat : 180 - dlat / 2;
lat0 = 90 - dlat / 2 : -dlat : -90 + dlat / 2;
sha_glq = @(xrho, xNq) tess2Vnm_glq(lon0, lat0, R1, R2, xrho, nmax, "a", a, "M", M, "Nq", xNq);
sha_ml = @(xrho) tess2Vnm_multLayer(lon0, lat0, R1, R2, xrho, nmax, "a", a, "M", M);
% fwd spe
tic
vnm_exp = tess2Vnm_exp(lon0, lat0, R1, R2, reshape(rho_exp, 1, 1, []), nmax, "a", a, "M", M, r0 = r0);
timeexp = toc;
re_exp = RE(vnm_exp(1, 3), vnm_exp_ref);
tic
vnm_poly = tess2Vnm_poly(lon0, lat0, R1, R2, reshape(rho_poly, 1, 1, []), nmax, "a", a, "M", M, r0 = r0);
timepoly = toc;
re_poly = RE(vnm_poly(1, 3), vnm_poly_ref);
tic
vnm_par = tess2Vnm_para(lon0, lat0, R1, R2, reshape(rho_par, 1, 1, []), nmax, "a", a, "M", M, r0 = r0);
timepar = toc;
re_par = RE(vnm_par(1, 3), vnm_par_ref);

tic
vnm_exp_ml = sha_ml(frho_exp);
timeexp_ml = toc;
re_exp_ml = RE(vnm_exp_ml(1, 3), vnm_exp_ref);
tic
vnm_poly_ml = sha_ml(frho_poly);
timepoly_ml = toc;
re_poly_ml = RE(vnm_poly_ml(1, 3), vnm_poly_ref);
tic
vnm_par_ml = sha_ml(frho_par);
timepar_ml = toc;
re_par_ml = RE(vnm_par_ml(1, 3), vnm_par_ref);
%%
Nq = 2 : 1 : 32; Numq = length(Nq);
re_poly_qlg = zeros(Numq, 1); Time_poly = re_poly_qlg;
for i = 1 : Numq
    tic
    vnm_glq = sha_glq(frho_poly, Nq(i));
    Time_poly(i) = toc;
    re_poly_qlg(i) = RE(vnm_glq(1, 3), vnm_poly_ref);
end

re_exp_qlg = zeros(Numq, 1); Time_exp = re_exp_qlg;
for i = 1 : Numq
    tic
    vnm_glq = sha_glq(frho_exp, Nq(i));
    Time_exp(i) = toc;
    re_exp_qlg(i) = RE(vnm_glq(1, 3), vnm_exp_ref);
end

re_par_qlg = zeros(Numq, 1); Time_par = re_par_qlg;
for i = 1 : Numq
    tic
    vnm_glq = sha_glq(frho_par, Nq(i));
    Time_par(i) = toc;
    re_par_qlg(i) = RE(vnm_glq(1, 3), vnm_par_ref);
end
save('shell_vd', 'Time_poly', 're_poly_qlg', 'Time_exp', 're_exp_qlg', 'Time_par', 're_par_qlg');

vnm_glq_poly = sha_glq(frho_poly, 5);
vnm_glq_exp = sha_glq(frho_exp, 11);
vnm_glq_par = sha_glq(frho_par, 20);
load('shell_vd', 'Time_poly', 're_poly_qlg', 'Time_exp', 're_exp_qlg', 'Time_par', 're_par_qlg');
%%
vnm = vnm_poly;
vnm(1, 3) = vnm(1, 3) - vnm_poly_ref;
dsn_poly = sqrtSigma(vnm);

vnm = vnm_poly_ml;
vnm(1, 3) = vnm(1, 3) - vnm_poly_ref;
dsn_ml_poly = sqrtSigma(vnm);

vnm = vnm_glq_poly;
vnm(1, 3) = vnm(1, 3) - vnm_poly_ref;
dsn_glq_poly = sqrtSigma(vnm);

vnm = vnm_exp;
vnm(1, 3) = vnm(1, 3) - vnm_exp_ref;
dsn_exp = sqrtSigma(vnm);

vnm = vnm_exp_ml;
vnm(1, 3) = vnm(1, 3) - vnm_exp_ref;
dsn_ml_exp = sqrtSigma(vnm);

vnm = vnm_glq_exp;
vnm(1, 3) = vnm(1, 3) - vnm_exp_ref;
dsn_glq_exp = sqrtSigma(vnm);

vnm = vnm_par;
vnm(1, 3) = vnm(1, 3) - vnm_par_ref;
dsn_par = sqrtSigma(vnm);

vnm = vnm_par_ml;
vnm(1, 3) = vnm(1, 3) - vnm_par_ref;
dsn_ml_par = sqrtSigma(vnm);

vnm = vnm_glq_par;
vnm(1, 3) = vnm(1, 3) - vnm_par_ref;
dsn_glq_par = sqrtSigma(vnm);

x = 0 : 1 : nmax;
%% plot
fontsize = 9;
linewidth = 1;
fig = figure;
outer = tiledlayout(fig, 30, 1, "Padding", 'tight', 'TileSpacing', 'tight');
inner = tiledlayout(3, 2, "Padding", 'tight', 'TileSpacing', 'tight', Parent=outer);
inner.Layout.TileSpan = [29, 1];

nexttile(inner)
plot(Nq, re_poly_qlg, "LineWidth", linewidth, "Color", [0 0.4470 0.7410])
hold on
line([2, 32], [re_poly, re_poly], 'LineStyle', '--', "LineWidth", linewidth, 'Color', [0.8500 0.3250 0.0980]);
line([2, 32], [re_poly_ml, re_poly_ml], 'LineStyle', '-.', "LineWidth", linewidth, 'Color', [0.4940 0.1840 0.5560]);
scatter(Nq(5 - 1) ,re_poly_qlg(5 - 1), 80, 'x', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
xlim([2, 32]);
xticks(2 : 6 : 32)
yticks(-16 : 4 : 0);
ylim([-16, 0]);
set(gca, 'YMinorTick', 'on');
ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex');
title('(a)', 'FontSize', fontsize, 'Interpreter', 'latex');
grid on
set(gca,'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
plot(Nq, Time_poly, "LineWidth", linewidth);
hold on
line([2, 32], [timepoly, timepoly], 'LineStyle', '--', "LineWidth", linewidth);
line([2, 32], [timepoly_ml, timepoly_ml], 'LineStyle', '-.', "LineWidth", linewidth, 'Color', [0.4940 0.1840 0.5560]);
scatter(Nq(5 - 1) , Time_poly(5 - 1), 80, 'x', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
xlim([2, 32]);
ylim([0, 40]);
yticks(0: 10: 40);
xticks(2 : 6 : 32)
set(gca, 'YMinorTick', 'on');
ylabel('Computation times (s)', 'Interpreter', 'latex');
title('(d)', 'FontSize', fontsize, 'Interpreter', 'latex');
grid on
set(gca,'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
plot(Nq, re_exp_qlg, "LineWidth", linewidth)
hold on
line([2, 32], [re_exp, re_exp], 'LineStyle', '--', "LineWidth", linewidth);
line([2, 32], [re_exp_ml, re_exp_ml], 'LineStyle', '-.', "LineWidth", linewidth, 'Color', [0.4940 0.1840 0.5560]);
scatter(Nq(11 - 1) , re_exp_qlg(11 - 1), 80, 'x', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
yticks(-16:4:0);
ylim([-16, 0]);
ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex');
xlim([2, 32]);
xticks(2 : 6 : 32)
yticks(-16 : 4 : 0);
ylim([-16, 0]);
set(gca, 'YMinorTick', 'on');
title('(b)', 'FontSize', fontsize, 'Interpreter', 'latex');
grid on
set(gca,'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
plot(Nq, Time_exp, "LineWidth", linewidth);
hold on
line([2, 32], [timeexp, timeexp], 'LineStyle', '--', "LineWidth", linewidth);
line([2, 32], [timeexp_ml, timeexp_ml], 'LineStyle', '-.', "LineWidth", linewidth, 'Color', [0.4940 0.1840 0.5560]);
scatter(Nq(11 - 1) , Time_exp(11 - 1), 80, 'x', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
xlim([2, 32]);
ylim([0, 40]);
yticks(0: 10: 40);
xticks(2 : 6 : 32)
title('(e)', 'FontSize', fontsize, 'Interpreter', 'latex');
grid on
ylabel('Computation times (s)', 'Interpreter', 'latex');
set(gca,'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
plot(Nq, re_par_qlg, "LineWidth", linewidth)
hold on
line([2, 32], [re_par, re_par], 'LineStyle', '--', "LineWidth", linewidth);
line([2, 32], [re_par_ml, re_par_ml], 'LineStyle', '-.', "LineWidth", linewidth, 'Color', [0.4940 0.1840 0.5560]);
scatter(Nq(20 - 1) , re_par_qlg(20 - 1), 80, 'x', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
yticks(-16:4:0);
ylim([-16, 0]);
xlim([2, 32]);
xticks(2 : 6 : 32)
yticks(-16 : 4 : 0);
ylim([-16, 0]);
set(gca, 'YMinorTick', 'on');
ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex');
title('(c)', 'FontSize', fontsize, 'Interpreter', 'latex');
grid on
set(gca,'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
plot(Nq, Time_par, "LineWidth", linewidth);
hold on
line([2, 32], [timepar, timepar], 'LineStyle', '--', "LineWidth", linewidth);
line([2, 32], [timepar_ml, timepar_ml], 'LineStyle', '-.', "LineWidth", linewidth, 'Color', [0.4940 0.1840 0.5560]);
scatter(Nq(20 - 1) , Time_par(20 - 1), 80, 'x', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
xlim([2, 32]);
ylim([0, 40]);
yticks(0: 10: 40);
xticks(2 : 6 : 32)
title('(f)', 'FontSize', fontsize, 'Interpreter', 'latex');
ylabel('Computation times (s)', 'Interpreter', 'latex');
xlabel(inner, 'The order of the quadrature $N_{q}$', 'Interpreter', 'latex', 'fontsize', fontsize);
grid on
set(gca,'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(outer)
plot(nan, nan, "LineWidth", linewidth);
hold on
line(nan, nan, 'LineStyle', '--', "LineWidth", linewidth);
line(nan, nan, 'LineStyle', '-.', "LineWidth", linewidth, 'Color', [0.4940 0.1840 0.5560]);

axis off;
scatter(nan , nan, 80, 'x', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
leg = legend('Numerical quadrature', 'Analytical integration', 'Block-wise constant', 'Optimal value of $N_q$', 'fontsize', fontsize, 'NumColumns', 2, 'Interpreter', 'latex');
leg.Layout.Tile = 'South';
leg.Box = "off";
leg.ItemTokenSize = [20, 20];
set(gcf,'Units','centimeters','Position', [1.3 3 13, 15]);
exportgraphics(gcf, '../manuscript/marked/shell_vd.pdf');
%%
sha_glq = @(xlon0, xlat0, xdlat, xdlon) tess2Vnm_glq(xlon0, xlat0, R1, R2, frho_poly, nmax, "a", a, "M", M, "Nq", [], "szm", [xdlat, xdlon]);
sha_ml = @(xlon0, xlat0, xdlat, xdlon) tess2Vnm_multLayer(xlon0, xlat0, R1, R2, frho_poly, nmax, "a", a, "M", M, "szm", [xdlat, xdlon]);
sha_poly = @(xlon0, xlat0, xdlat, xdlon) tess2Vnm_poly(xlon0, xlat0, R1, R2, reshape(rho_poly, 1, 1, []), nmax, "a", a, "M", M, 'r0', r0, 'szm', [xdlat, xdlon]);
dlat = 30; dlon = dlat;
lat0 = 90 - dlat / 2 : -dlat : -90 + dlat / 2;
lon0 = -180 + dlon / 2 : dlon: 180 - dlon / 2;
tic
vnm_30_as = sha_poly(lon0, lat0, dlat, dlon);
time_30_as = toc;
tic
vnm_30_glq = sha_glq(lon0, lat0, dlat, dlon);
time_30_glq = toc;
tic
vnm_30_ml = sha_ml(lon0, lat0, dlat, dlon);
time_30_ml = toc;
re_30_as = RE(vnm_30_as(1, 3), vnm_poly_ref);
re_30_glq = RE(vnm_30_glq(1, 3), vnm_poly_ref);
re_30_ml = RE(vnm_30_ml(1, 3), vnm_poly_ref);

dlat = 10; dlon = dlat;
lat0 = 90 - dlat / 2 : -dlat : -90 + dlat / 2;
lon0 = -180 + dlon / 2 : dlon: 180 - dlon / 2;
tic
vnm_10_as = sha_poly(lon0, lat0, dlat, dlon);
time_10_as = toc;
tic
vnm_10_glq = sha_glq(lon0, lat0, dlat, dlon);
time_10_glq = toc;
tic
vnm_10_ml = sha_ml(lon0, lat0, dlat, dlon);
time_10_ml = toc;
re_10_as = RE(vnm_10_as(1, 3), vnm_poly_ref);
re_10_glq = RE(vnm_10_glq(1, 3), vnm_poly_ref);
re_10_ml = RE(vnm_10_ml(1, 3), vnm_poly_ref);

dlat = 1; dlon = dlat;
lat0 = 90 - dlat / 2 : -dlat : -90 + dlat / 2;
lon0 = -180 + dlon / 2 : dlon: 180 - dlon / 2;
tic
vnm_1_as = sha_poly(lon0, lat0, dlat, dlon);
time_1_as = toc;
tic
vnm_1_glq = sha_glq(lon0, lat0, dlat, dlon);
time_1_glq = toc;
tic
vnm_1_ml = sha_ml(lon0, lat0, dlat, dlon);
time_1_ml = toc;
re_1_as = RE(vnm_1_as(1, 3), vnm_poly_ref);
re_1_glq = RE(vnm_1_glq(1, 3), vnm_poly_ref);
re_1_ml = RE(vnm_1_ml(1, 3), vnm_poly_ref);

dlat = 1 / 2; dlon = dlat;
lat0 = 90 - dlat / 2 : -dlat : -90 + dlat / 2;
lon0 = -180 + dlon / 2 : dlon: 180 - dlon / 2;
tic
vnm_1_2_as = sha_poly(lon0, lat0, dlat, dlon);
time_1_2_as = toc;
tic
vnm_1_2_glq = sha_glq(lon0, lat0, dlat, dlon);
time_1_2_glq = toc;
tic
vnm_1_2_ml = sha_ml(lon0, lat0, dlat, dlon);
time_1_2_ml = toc;
re_1_2_as = RE(vnm_1_2_as(1, 3), vnm_poly_ref);
re_1_2_glq = RE(vnm_1_2_glq(1, 3), vnm_poly_ref);
re_1_2_ml = RE(vnm_1_2_ml(1, 3), vnm_poly_ref);

dlat = 1 / 12; dlon = dlat;
lat0 = 90 - dlat / 2 : -dlat : -90 + dlat / 2;
lon0 = -180 + dlon / 2 : dlon: 180 - dlon / 2;
tic
vnm_1_12_as = sha_poly(lon0, lat0, dlat, dlon);
time_1_12_as = toc;
tic
vnm_1_12_glq = sha_glq(lon0, lat0, dlat, dlon);
time_1_12_glq = toc;
tic
vnm_1_12_ml = sha_ml(lon0, lat0, dlat, dlon);
time_1_12_ml = toc;

re_1_12_as = RE(vnm_1_12_as(1, 3), vnm_poly_ref);
re_1_12_glq = RE(vnm_1_12_glq(1, 3), vnm_poly_ref);
re_1_12_ml = RE(vnm_1_12_ml(1, 3), vnm_poly_ref);