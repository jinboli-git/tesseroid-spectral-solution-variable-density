clear
close all
clc
% parameter
R2 = 6371e3;
R1 = 6361e3;
h = 10 * 1e3;
r = R2 + h;
G = 6.672e-11;
M = 5.965e24;
a = R2;
rho0 = 2670;
clon = 0; clat = 0;
fref = @(xh) shell_geff_cons(R1, R2, R2 + xh, rho0, G);
RE = @(x, xref) log10(abs(x ./ xref - 1));
AE = @(x, xref) log10(abs(x - xref));
vnm_ref = shell_C00_cons(R1, R2, rho0, M);
sqrtSigma = @(x) log10( sqrt(degVar(x)) );
shs = @(x, xnmax, xh, xtype) shs_grid(x, clon, clat, xnmax, 0, "height", xh, "fldTyp", xtype, "GM", G * M, "a", a);
shs2 = @(x, xclat, xnmax, xh, xtype) shs_grid(x, clon, xclat, xnmax, 0, "height", xh, "fldTyp", xtype, "GM", G * M, "a", a);
sha = @(xlon0, xlat0, xnmax, xdlat, xdlon) tess2Vnm(xlon0, xlat0, R1, R2, rho0, xnmax, "a", a, "M", M, "szm", [xdlat, xdlon]);
sha_fft = @(xlon0, xlat0, xnmax, xdlat, xdlon) tess2Vnm_fft(xlon0, xlat0, R1, R2, rho0, xnmax, "a", a, "M", M, "szm", [xdlat, xdlon]);
nmax = 2190;
%% dsigman
rng("default");
sd = rng;
dlat = 30; dlon = dlat;
lat0 = 90 - dlat / 2 : -dlat : -90 + dlat / 2;
lon0 = -180 + dlon / 2 : dlon: 180 - dlon / 2;
lon0 = lon0(randperm(length(lon0)));
vnm_30 = sha(lon0, lat0, nmax, dlat, dlon);

dlat = 10; dlon = dlat;
lat0 = 90 - dlat / 2 : -dlat : -90 + dlat / 2;
lon0 = -180 + dlon / 2 : dlon: 180 - dlon / 2;
lon0 = lon0(randperm(length(lon0)));
vnm_10 = sha(lon0, lat0, nmax, dlat, dlon);

dlat = 1; dlon = dlat;
lat0 = 90 - dlat / 2 : -dlat : -90 + dlat / 2;
lon0 = -180 + dlon / 2 : dlon: 180 - dlon / 2;
lon0 = lon0(randperm(length(lon0)));
vnm_1 = sha(lon0, lat0, nmax, dlat, dlon);

dlat = 1 / 2; dlon = dlat;
lat0 = 90 - dlat / 2 : -dlat : -90 + dlat / 2;
lon0 = -180 + dlon / 2 : dlon: 180 - dlon / 2;
lon0 = lon0(randperm(length(lon0)));
vnm_1_2 = sha(lon0, lat0, nmax, dlat, dlon);
% 
dlat = 1 / 12; dlon = dlat;
lat0 = 90 - dlat / 2 : -dlat : -90 + dlat / 2;
lon0 = -180 + dlon / 2 : dlon: 180 - dlon / 2;
lon0 = lon0(randperm(length(lon0)));
vnm_1_12 = sha_fft(lon0, lat0, nmax, dlat, dlon);

% save('test_shell_dsn.mat', 'vnm_30', 'vnm_10', 'vnm_1', 'vnm_1_2', 'vnm_1_12');
%%
% load 'test_shell_dsn.mat'
vnm = vnm_30;
re_c00_30 = RE(vnm(1, 3), vnm_ref); ae_c00_30 = AE(vnm(1, 3), vnm_ref);
vnm(1, 3) = vnm(1, 3) - vnm_ref;
dsn_30 = sqrtSigma(vnm);

vnm = vnm_10;
re_c00_10 = RE(vnm(1, 3), vnm_ref); ae_c00_10 = AE(vnm(1, 3), vnm_ref);
vnm(1, 3) = vnm(1, 3) - vnm_ref;
dsn_10 = sqrtSigma(vnm);

vnm = vnm_1;
re_c00_1 = RE(vnm(1, 3), vnm_ref); ae_c00_1 = AE(vnm(1, 3), vnm_ref);
vnm(1, 3) = vnm(1, 3) - vnm_ref;
dsn_1 = sqrtSigma(vnm);

vnm = vnm_1_2;
re_c00_1_2 = RE(vnm(1, 3), vnm_ref); ae_c00_1_2 = AE(vnm(1, 3), vnm_ref);
vnm(1, 3) = vnm(1, 3) - vnm_ref;
dsn_1_2 = sqrtSigma(vnm);

vnm = vnm_1_12;
re_c00_1_12 = RE(vnm(1, 3), vnm_ref); ae_c00_1_12 = AE(vnm(1, 3), vnm_ref);
vnm(1, 3) = vnm(1, 3) - vnm_ref;
dsn_1_12 = sqrtSigma(vnm);

x = 0 : 1 : nmax;
dsn_30 = dsn_30(x + 1);
dsn_10 = dsn_10(x + 1);
dsn_1 = dsn_1(x + 1);
dsn_1_2 = dsn_1_2(x + 1);
dsn_1_12 = dsn_1_12(x + 1);
%% plot
fontsize = 9;
cmp = [0.85, 0.10, 0.10; 0.00, 0.45, 0.74; 0.20, 0.60, 0.20; 1.00, 0.60, 0.00; 0.50, 0.20, 0.70; [255 215 0]/255; 0.55, 0.27, 0.07; 1.00, 0.60, 0.80; 0.6, 0.6, 0.6];
fig = figure;
linewidth = 0.5;
tiledlayout(1, 1, "Padding", "compact", 'TileSpacing', 'compact', Parent=fig);
nexttile
plot(x, dsn_30, 'LineWidth', linewidth, 'Color', cmp(1, :));
hold on
plot(x, dsn_10, 'LineWidth', linewidth, 'Color', cmp(2, :));
plot(x, dsn_1, 'LineWidth', linewidth, 'Color', cmp(3, :));
plot(x, dsn_1_2, 'LineWidth', linewidth, 'Color', cmp(4, :));
plot(x, dsn_1_12, 'LineWidth', linewidth, 'Color', cmp(5, :));

set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);
xlim([0, 2200]);
ylim([-28, -16]);
xticks(0:200:2200);
yticks(-28:2:-16);
grid on
ylabel('Square root of error degree variance in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize);
xlabel('Spherical harmonic d/o', 'Interpreter', 'latex', 'fontsize', fontsize);

leg = legend('$\Delta\theta^\prime=\Delta\lambda^\prime=30^\circ$', '$10^\circ$', '$1^\circ$', '$1/2^\circ$', '$5^\prime$', 'fontsize', fontsize, 'NumColumns', 6, 'Interpreter', 'latex');
leg.Layout.Tile = 'South';
leg.Box = "off";

set(gcf,'Units','centimeters','Position',[1.3, 3, 13, 8]);
exportgraphics(gcf, '../manuscript/JG/figures/shell_SEDV.pdf');
%%
vnm = vnm_1_12;
H = 0 : 10 : 500;
H = H * 1e3; H(1) = 10;
num = length(H);

ae_v = zeros(num, 1);
ae_vz = zeros(num, 1);
ae_vzz = zeros(num, 1);
ae_vxx = zeros(num, 1);
ae_vyy = zeros(num, 1);
ae_vzzz = zeros(num, 1);
ae_vxxz = zeros(num, 1);
ae_vyyz = zeros(num, 1);
L1 = zeros(num, 1);
L2 = zeros(num, 1);

re_v = zeros(num, 1);
re_vz = zeros(num, 1);
re_vzz = zeros(num, 1);
re_vxx = zeros(num, 1);
re_vyy = zeros(num, 1);
re_vzzz = zeros(num, 1);
re_vxxz = zeros(num, 1);
re_vyyz = zeros(num, 1);

for i = 1 : num
    h = H(i);
    [v_ref, vz_ref, vzz_ref, vxx_ref, vzzz_ref, vxxz_ref] = fref(h);

    v_cal = shs(vnm, nmax, h, 'v');
    vz_cal = shs(vnm, nmax, h, 'vz') / 1e5;
    vzz_cal = shs(vnm, nmax, h, 'vzz') / 1e9;
    vxx_cal = shs(vnm, nmax, h, 'vxx') / 1e9;
    vyy_cal = shs(vnm, nmax, h, 'vyy') / 1e9;
    vzzz_cal = shs(vnm, nmax, h, 'vzzz') / 1e15;
    vxxz_cal = shs(vnm, nmax, h, 'vxxz') / 1e15;
    vyyz_cal = shs(vnm, nmax, h, 'vyyz') / 1e15;

    ae_v(i) = AE(v_cal, v_ref);
    ae_vz(i) = AE(vz_cal, vz_ref);
    ae_vzz(i) = AE(vzz_cal, vzz_ref);
    ae_vxx(i) = AE(vxx_cal, vxx_ref);
    ae_vyy(i) = AE(vyy_cal, vxx_ref);
    ae_vzzz(i) = AE(vzzz_cal, vzzz_ref);
    ae_vxxz(i) = AE(vxxz_cal, vxxz_ref);
    ae_vyyz(i) = AE(vyyz_cal, vxxz_ref);
    L1(i) = log10(abs(vxx_cal + vyy_cal + vzz_cal));
    L2(i) = log10(abs(vzzz_cal + vxxz_cal + vyyz_cal));

    re_v(i) = RE(v_cal, v_ref);
    re_vz(i) = RE(vz_cal, vz_ref);
    re_vzz(i) = RE(vzz_cal, vzz_ref);
    re_vxx(i) = RE(vxx_cal, vxx_ref);
    re_vyy(i) = RE(vyy_cal, vxx_ref);
    re_vzzz(i) = RE(vzzz_cal, vzzz_ref);
    re_vxxz(i) = RE(vxxz_cal, vxxz_ref);
    re_vyyz(i) = RE(vyyz_cal, vxxz_ref);
end
H = H' / 1e3;
%% plot
fontsize = 9;
fig = figure;
outer = tiledlayout(fig, 15, 1, "Padding", 'tight', 'TileSpacing', 'tight');
inner = tiledlayout(1, 2, "Padding", 'tight', 'TileSpacing', 'tight', Parent=outer);
inner.Layout.TileSpan = [14, 1];
mksz = 24;
linewidth = 0.8;
nexttile(inner)
scatter(H, ae_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(H, ae_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(H, ae_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(H, ae_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(H, ae_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(H, ae_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(H, ae_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(H, ae_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
scatter(H, L1, mksz, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
scatter(H, L2, mksz, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
xlim([0, 500])
grid on
box on
title('(a)');
ylabel('Absolute errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
scatter(H, re_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(H, re_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(H, re_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(H, re_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(H, re_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(H, re_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(H, re_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(H, re_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);

ylim([-15, -7]);
xlim([0, 500])
yticks(-15:2:-7)
ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
grid on
box on
title('(b)');
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);
xlabel(inner, 'Observation height $h$ (km)', 'Interpreter', 'latex', 'fontsize', fontsize)

nexttile(outer)

scatter(nan, nan, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(nan, nan, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
scatter(nan, nan, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
scatter(nan, nan, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);

leg = legend('$V$', '$V_z$', '$V_{zz}$', '$V_{xx}$', '$V_{yy}$', '$V_{zzz}$', '$V_{xxz}$', '$V_{yyz}$', '$\delta L_1$','$\delta L_2$', 'fontsize', fontsize, 'NumColumns', 6, 'Interpreter', 'latex');
leg.Box = "off";
leg.Layout.Tile = 'South';

axis off;

set(gcf,'Units','centimeters','Position', [1.3 3 13, 8]);
exportgraphics(gcf, '../manuscript/JG/figures/shell_H.pdf');
%%
h = 10 * 1e3;
[v_ref, vz_ref, vzz_ref, vxx_ref, vzzz_ref, vxxz_ref] = fref(h);
Clat1 = 70 : 1/12 : 90;
v_cal = shs2(vnm, Clat1, nmax, h, 'v');
vz_cal = shs2(vnm, Clat1, nmax, h, 'vz') / 1e5;
vzz_cal = shs2(vnm, Clat1, nmax, h, 'vzz') / 1e9;
vxx_cal = shs2(vnm, Clat1, nmax, h, 'vxx') / 1e9;
vyy_cal = shs2(vnm, Clat1, nmax, h, 'vyy') / 1e9;
vzzz_cal = shs2(vnm, Clat1, nmax, h, 'vzzz') / 1e15;
vxxz_cal = shs2(vnm, Clat1, nmax, h, 'vxxz') / 1e15;
vyyz_cal = shs2(vnm, Clat1, nmax, h, 'vyyz') / 1e15;

ae_clat1_v = AE(v_cal, v_ref);
ae_clat1_vz = AE(vz_cal, vz_ref);
ae_clat1_vzz = AE(vzz_cal, vzz_ref);
ae_clat1_vxx = AE(vxx_cal, vxx_ref);
ae_clat1_vyy = AE(vyy_cal, vxx_ref);
ae_clat1_vzzz = AE(vzzz_cal, vzzz_ref);
ae_clat1_vxxz = AE(vxxz_cal, vxxz_ref);
ae_clat1_vyyz = AE(vyyz_cal, vxxz_ref);
L1_clat1 = log10(abs(vxx_cal + vyy_cal + vzz_cal));
L2_clat1 = log10(abs(vzzz_cal + vxxz_cal + vyyz_cal));

re_clat1_v = RE(v_cal, v_ref);
re_clat1_vz = RE(vz_cal, vz_ref);
re_clat1_vzz = RE(vzz_cal, vzz_ref);
re_clat1_vxx = RE(vxx_cal, vxx_ref);
re_clat1_vyy = RE(vyy_cal, vxx_ref);
re_clat1_vzzz = RE(vzzz_cal, vzzz_ref);
re_clat1_vxxz = RE(vxxz_cal, vxxz_ref);
re_clat1_vyyz = RE(vyyz_cal, vxxz_ref);

Clat2 = -70 : -1/12 : -90;
v_cal = shs2(vnm, Clat2, nmax, h, 'v');
vz_cal = shs2(vnm, Clat2, nmax, h, 'vz') / 1e5;
vzz_cal = shs2(vnm, Clat2, nmax, h, 'vzz') / 1e9;
vxx_cal = shs2(vnm, Clat2, nmax, h, 'vxx') / 1e9;
vyy_cal = shs2(vnm, Clat2, nmax, h, 'vyy') / 1e9;
vzzz_cal = shs2(vnm, Clat2, nmax, h, 'vzzz') / 1e15;
vxxz_cal = shs2(vnm, Clat2, nmax, h, 'vxxz') / 1e15;
vyyz_cal = shs2(vnm, Clat2, nmax, h, 'vyyz') / 1e15;

ae_clat2_v = AE(v_cal, v_ref);
ae_clat2_vz = AE(vz_cal, vz_ref);
ae_clat2_vzz = AE(vzz_cal, vzz_ref);
ae_clat2_vxx = AE(vxx_cal, vxx_ref);
ae_clat2_vyy = AE(vyy_cal, vxx_ref);
ae_clat2_vzzz = AE(vzzz_cal, vzzz_ref);
ae_clat2_vxxz = AE(vxxz_cal, vxxz_ref);
ae_clat2_vyyz = AE(vyyz_cal, vxxz_ref);
L1_clat2 = log10(abs(vxx_cal + vyy_cal + vzz_cal));
L2_clat2 = log10(abs(vzzz_cal + vxxz_cal + vyyz_cal));

re_clat2_v = RE(v_cal, v_ref);
re_clat2_vz = RE(vz_cal, vz_ref);
re_clat2_vzz = RE(vzz_cal, vzz_ref);
re_clat2_vxx = RE(vxx_cal, vxx_ref);
re_clat2_vyy = RE(vyy_cal, vxx_ref);
re_clat2_vzzz = RE(vzzz_cal, vzzz_ref);
re_clat2_vxxz = RE(vxxz_cal, vxxz_ref);
re_clat2_vyyz = RE(vyyz_cal, vxxz_ref);
%% plot
Clat1 = 90 - Clat1;
Clat2 = 90 - Clat2;
fontsize = 9;
fig = figure;
outer = tiledlayout(fig, 30, 1, "Padding", 'tight', 'TileSpacing', 'tight');
inner = tiledlayout(2, 2, 'Padding', "compact", 'TileSpacing', 'compact', Parent=outer);
inner.Layout.TileSpan = [29, 1];

mksz = 8;
linewidth = 0.8;

nexttile(inner)
scatter(Clat1, ae_clat1_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(Clat1, ae_clat1_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(Clat1, ae_clat1_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(Clat1, ae_clat1_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(Clat1, ae_clat1_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(Clat1, ae_clat1_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(Clat1, ae_clat1_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(Clat1, ae_clat1_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
scatter(Clat1, L1_clat1, mksz, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
scatter(Clat1, L2_clat1, mksz, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
xlim([0, 20]);
xticks(0 : 5 : 20);
ylim([-35, -5]);
grid on
box on
title('(a)');
ylabel('Absolute errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize);
set(gca, 'XDir', 'reverse');
set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabel', {'$0^\circ$', '$5^\circ$', '$10^\circ$', '$15^\circ$', '$20^\circ$'});
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
scatter(Clat1, re_clat1_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(Clat1, re_clat1_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(Clat1, re_clat1_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(Clat1, re_clat1_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(Clat1, re_clat1_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(Clat1, re_clat1_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(Clat1, re_clat1_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(Clat1, re_clat1_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
xlim([0, 20]);
xticks(0 : 5 : 20);
ylim([-16, -6]);
yticks(-16 : 2 : -6)
ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
grid on
box on
title('(b)');
set(gca, 'XDir', 'reverse');
set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabel', {'$0^\circ$', '$5^\circ$', '$10^\circ$', '$15^\circ$', '$20^\circ$'});
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);
xlabel(inner, 'Computation point colatitude $\theta$', 'Interpreter', 'latex', 'fontsize', fontsize)

nexttile(inner)
scatter(Clat2, ae_clat2_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(Clat2, ae_clat2_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(Clat2, ae_clat2_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(Clat2, ae_clat2_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(Clat2, ae_clat2_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(Clat2, ae_clat2_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(Clat2, ae_clat2_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(Clat2, ae_clat2_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
scatter(Clat2, L1_clat2, mksz, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
scatter(Clat2, L2_clat2, mksz, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
xlim([160, 180]);
xticks(160 : 5 : 180);
ylim([-35, -5]);
grid on
box on
title('(c)');
ylabel('Absolute errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize);
set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabel', {'$160^\circ$', '$165^\circ$', '$170^\circ$', '$175^\circ$', '$180^\circ$'});
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
scatter(Clat2, re_clat2_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(Clat2, re_clat2_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(Clat2, re_clat2_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(Clat2, re_clat2_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(Clat2, re_clat2_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(Clat2, re_clat2_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(Clat2, re_clat2_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(Clat2, re_clat2_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
xlim([160, 180]);
xticks(160 : 5 : 180);
ylim([-16, -6]);
yticks(-16 : 2 : -6)
ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
grid on
box on
title('(d)');
set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabel', {'$160^\circ$', '$165^\circ$', '$170^\circ$', '$175^\circ$', '$180^\circ$'});
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(outer)
scatter(nan, nan, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(nan, nan, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
scatter(nan, nan, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
scatter(nan, nan, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);

leg = legend('$V$', '$V_z$', '$V_{zz}$', '$V_{xx}$', '$V_{yy}$', '$V_{zzz}$', '$V_{xxz}$', '$V_{yyz}$', '$\delta L_1$','$\delta L_2$', 'fontsize', fontsize, 'NumColumns', 6, 'Interpreter', 'latex');
leg.Box = "off";
leg.Layout.Tile = 'South';
axis off;

set(gcf,'Units','centimeters','Position', [1.3 3 13, 13]);
exportgraphics(gcf, '../manuscript/JG/figures/shell_Clat.pdf');