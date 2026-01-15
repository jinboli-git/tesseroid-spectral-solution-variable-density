clear
close all
clc
% parameter
R2 = 6371e3;
R1 = 6361e3;
G = 6.672e-11;
M = 5.965e24;
a = R2;
rho = 2670;

lon0 = 0;
clon = 0;
clat = 90;
RE = @(x, xref) log10(abs(x / xref - 1));
AE = @(x, xref) log10(abs(x - xref));
sha = @(xlon0, xlat0, xnmax, xdlat, xdlon) tess2Vnm(xlon0, xlat0, R1, R2, rho, xnmax, "a", a, "M", M, szm = [xdlat, xdlon]);
shs = @(x, xnmax, xh, xtyp) shs_grid(x, clon, clat, xnmax, 0, "height", xh, "fldTyp", xtyp, "GM", G * M, "a", a);
dsigma = @(x) degVar(x);
pow = @(x, xh, xtype) power_log10(x, G * M, R2 + xh, a, xtype);
%% d/o vs RVD: sz
lat0 = 60; 
dlat = 30;
nmax = 2190;
vnm_30 = sha(lon0, lat0, nmax, dlat, dlat);
dlat = 10;
vnm_10 = sha(lon0, lat0, nmax, dlat, dlat);
dlat = 1;
vnm_1 = sha(lon0, lat0, nmax, dlat, dlat);
dlat = 1 / 2;
vnm_1_2 = sha(lon0, lat0, nmax, dlat, dlat);
dlat = 1 / 12;
vnm_1_12 = sha(lon0, lat0, nmax, dlat, dlat);
dsigma = @(x) log10( dsigma(x) / sum(dsigma(x)) );
dsn_1_12 = dsigma(vnm_1_12);
dsn_1_2 = dsigma(vnm_1_2);
dsn_30 = dsigma(vnm_30);
dsn_1 = dsigma(vnm_1);
dsn_10 = dsigma(vnm_10);
%% plot
fontsize = 9;
cmp = [0.85, 0.10, 0.10; 0.00, 0.45, 0.74; 0.20, 0.60, 0.20; 1.00, 0.60, 0.00; 0.50, 0.20, 0.70; [255 215 0] / 255; 0.55, 0.27, 0.07; 1.00, 0.60, 0.80; 0.6, 0.6, 0.6; 0, 0, 0];
fig = figure;
tiledlayout(1, 1, "Padding", "compact", 'TileSpacing', 'compact', Parent=fig);
nexttile
x = 0 : 1 : nmax;
plot(x, dsn_30, 'LineWidth', 1, 'Color', cmp(1, :));
hold on
plot(x, dsn_10, 'LineWidth', 1, 'Color', cmp(2, :));
plot(x, dsn_1, 'LineWidth', 1, 'Color', cmp(3, :));
plot(x, dsn_1_2, 'LineWidth', 1, 'Color', cmp(4, :));
plot(x, dsn_1_12, 'LineWidth', 1, 'Color', cmp(5, :));

set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);
xlim([0, 2200]);
xticks(0 : 200 : 2200);
grid on
ylabel('Square root of error degree variance in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize);
xlabel('Spherical harmonic d/o', 'Interpreter', 'latex', 'fontsize', fontsize);

leg = legend('$\Delta\theta^\prime=\Delta\lambda^\prime=30^\circ$', '$10^\circ$', '$1^\circ$', '$1/2^\circ$', '$5^\prime$', 'fontsize', fontsize, 'NumColumns', 5, 'Interpreter', 'latex', 'Orientation', 'horizontal');
leg.Layout.Tile = 'South';
leg.Box = "off";

ylabel('Normalized degree variances in $\log_{10}$', 'Interpreter', 'latex');
set(gca, 'FontSize', fontsize);

set(gcf,'Units','centimeters','Position',[1.3, 3, 13, 8]);
exportgraphics(gcf, '../manuscript/marked/tess_NDV.pdf');
%%
dlat = 1;
vnm = sha(lon0, lat0, 10000, dlat, dlat);

hh = 260 * 1e3;
pv_1_260 = pow(vnm, hh, 'v');
pvz_1_260 = pow(vnm, hh, 'vz');
pvzz_1_260 = pow(vnm, hh, 'vzz');
pvzzz_1_260 = pow(vnm, hh, 'vzzz');

hh = 10 * 1e3;
pv_1_10000 = pow(vnm, hh, 'v');
pvz_1_10000 = pow(vnm, hh, 'vz');
pvzz_1_10000 = pow(vnm, hh, 'vzz');
pvzzz_1_10000 = pow(vnm, hh, 'vzzz');

dlat = 1/12;
vnm_1_2 = sha(lon0, lat0, 10000, dlat, dlat);

hh = 260 * 1e3;
pv_12_260 = pow(vnm_1_2, hh, 'v');
pvz_12_260 = pow(vnm_1_2, hh, 'vz');
pvzz_12_260 = pow(vnm_1_2, hh, 'vzz');
pvzzz_12_260 = pow(vnm_1_2, hh, 'vzzz');

hh = 10 * 1e3;
pv_12_10000 = pow(vnm_1_2, hh, 'v');
pvz_12_10000 = pow(vnm_1_2, hh, 'vz');
pvzz_12_10000 = pow(vnm_1_2, hh, 'vzz');
pvzzz_12_10000 = pow(vnm_1_2, hh, 'vzzz');
%% plot
fig = figure;
outer = tiledlayout(fig, 30, 1, "Padding", 'tight', 'TileSpacing', 'tight');
inner = tiledlayout(3, 1, "Padding", 'tight', 'TileSpacing', 'tight', Parent=outer);
inner.Layout.TileSpan = [29, 1];

x = 0 : 1 : 10000;
nexttile(inner, 1, [2, 1]);

plot(x, pv_1_260, 'LineWidth', 1, 'Color', cmp(1, :), 'LineStyle', '--');
hold on
plot(x, pvz_1_260, 'LineWidth', 1, 'Color', cmp(2, :), 'LineStyle', '--');
plot(x, pvzz_1_260, 'LineWidth', 1, 'Color', cmp(3, :), 'LineStyle', '--');
plot(x, pvzzz_1_260, 'LineWidth', 1, 'Color', cmp(5, :), 'LineStyle', '--');
plot(x, pv_12_260, 'LineWidth', 1, 'Color', cmp(1, :));
plot(x, pvz_12_260, 'LineWidth', 1, 'Color', cmp(2, :));
plot(x, pvzz_12_260, 'LineWidth', 1, 'Color', cmp(3, :));
plot(x, pvzzz_12_260, 'LineWidth', 1, 'Color', cmp(5, :));

set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);
ylim([-500, 0]);
xlim([0, 10000]);
xticks(0 : 2000 : 10000);
yticks(-500 : 100 : 0);
grid on
title('(a) $h=260\ \mathrm{km}$', 'Interpreter', 'latex', 'fontsize', fontsize);
ylabel(inner, 'Normalized power spectra in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize);
xlabel(inner, 'Spherical harmonic d/o', 'Interpreter', 'latex', 'fontsize', fontsize);
set(gca, 'FontSize', fontsize);
rectangle('Position', [0, -30, 600, 30], 'LineWidth', 1.2);
annotation('arrow', [0.14, 0.25], [0.92, 0.78], 'Color', [0.6, 0.6, 0.6, 0.6]);

ax_inset = axes('Units', 'centimeters', 'Position', [2.1, 7.3, 5, 4], Parent=fig);
box on
plot(x, pv_1_260, 'LineWidth', 1, 'Color', cmp(1, :), 'LineStyle', '--'); 
hold on
plot(x, pvz_1_260, 'LineWidth', 1, 'Color', cmp(2, :), 'LineStyle', '--');
plot(x, pvzz_1_260, 'LineWidth', 1, 'Color', cmp(3, :), 'LineStyle', '--');
plot(x, pvzzz_1_260, 'LineWidth', 1, 'Color', cmp(5, :), 'LineStyle', '--');
plot(x, pv_12_260, 'LineWidth', 1, 'Color', cmp(1, :));
plot(x, pvz_12_260, 'LineWidth', 1, 'Color', cmp(2, :));
plot(x, pvzz_12_260, 'LineWidth', 1, 'Color', cmp(3, :));
plot(x, pvzzz_12_260, 'LineWidth', 1, 'Color', cmp(5, :));

ylim([-30, 0]);
xlim([0, 600]);
xticks(0 : 200 : 600);
grid on
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

x = 0 : 1 : 10000;
nexttile(inner)
plot(x, pv_1_10000, 'LineWidth', 1, 'Color', cmp(1, :), 'LineStyle', '--');
hold on
plot(x, pvz_1_10000, 'LineWidth', 1, 'Color', cmp(2, :), 'LineStyle', '--');
plot(x, pvzz_1_10000, 'LineWidth', 1, 'Color', cmp(3, :), 'LineStyle', '--');
plot(x, pvzzz_1_10000, 'LineWidth', 1, 'Color', cmp(5, :), 'LineStyle', '--');
plot(x, pv_12_10000, 'LineWidth', 1, 'Color', cmp(1, :));
plot(x, pvz_12_10000, 'LineWidth', 1, 'Color', cmp(2, :));
plot(x, pvzz_12_10000, 'LineWidth', 1, 'Color', cmp(3, :));
plot(x, pvzzz_12_10000, 'LineWidth', 1, 'Color', cmp(5, :));
ylim([-30, 0]);
xlim([0, 10000]);
xticks(0 : 2000 : 10000);
grid on
title('(b) $h=10\ \mathrm{km}$', 'Interpreter', 'latex', 'fontsize', fontsize);
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(outer)
plot(nan, nan, 'LineWidth', 1, 'Color', cmp(1, :), 'LineStyle', '--');
hold on
plot(nan, nan, 'LineWidth', 1, 'Color', cmp(2, :), 'LineStyle', '--');
plot(nan, nan, 'LineWidth', 1, 'Color', cmp(3, :), 'LineStyle', '--');
plot(nan, nan, 'LineWidth', 1, 'Color', cmp(5, :), 'LineStyle', '--');
plot(nan, nan, 'LineWidth', 1, 'Color', cmp(1, :));
plot(nan, nan, 'LineWidth', 1, 'Color', cmp(2, :));
plot(nan, nan, 'LineWidth', 1, 'Color', cmp(3, :));
plot(nan, nan, 'LineWidth', 1, 'Color', cmp(5, :));
axis off;

leg = legend('$S/1^\circ$', '$S_{z}/1^\circ$', '$S_{zz}/1^\circ$', '$S_{zzz}/1^\circ$', '$S/5^\prime$', '$S_{z}/5^\prime$', '$S_{zz}/5^\prime$', '$S_{zzz}/5^\prime$', 'fontsize', fontsize, 'NumColumns', 4, 'Interpreter', 'latex', 'Orientation', 'horizontal');
leg.ItemTokenSize = [15, 15];
leg.Layout.Tile = 'South';
leg.Box = "off";

set(gcf,'Units','centimeters','Position',[1.3, 3, 13, 15]);
exportgraphics(gcf, '../manuscript/marked/tess_NPS.pdf', 'ContentType','vector');
%% h
dlat = 1;
nmax = 180;
vnm = sha(lon0, lat0, nmax, dlat, dlat);
H = 0 : 10 : 500;
H = H * 1e3; H(1) = 10;
num = length(H);

ae_v_1 = zeros(num, 1);
ae_vz_1 = zeros(num, 1);
ae_vzz_1 = zeros(num, 1);
ae_vxx_1 = zeros(num, 1);
ae_vyy_1 = zeros(num, 1);
ae_vzzz_1 = zeros(num, 1);
ae_vxxz_1 = zeros(num, 1);
ae_vyyz_1 = zeros(num, 1);
L1_1 = zeros(num, 1);
L2_1 = zeros(num, 1);

re_v_1 = zeros(num, 1);
re_vz_1 = zeros(num, 1);
re_vzz_1 = zeros(num, 1);
re_vxx_1 = zeros(num, 1);
re_vyy_1 = zeros(num, 1);
re_vzzz_1 = zeros(num, 1);
re_vxxz_1 = zeros(num, 1);
re_vyyz_1 = zeros(num, 1);

load ../data/vref_h_1.mat

for i = 1 : num
    h = H(i);
    v_cal = shs(vnm, nmax, h, 'v');
    vz_cal = shs(vnm, nmax, h, 'vz');
    vzz_cal = shs(vnm, nmax, h, 'vzz');
    vxx_cal = shs(vnm, nmax, h, 'vxx');
    vyy_cal = shs(vnm, nmax, h, 'vyy');
    vzzz_cal = shs(vnm, nmax, h, 'vzzz');
    vxxz_cal = shs(vnm, nmax, h, 'vxxz');
    vyyz_cal = shs(vnm, nmax, h, 'vyyz');

    ae_v_1(i) = AE(v_cal, v_ref(i));
    ae_vz_1(i) = AE(vz_cal, vz_ref(i));
    ae_vzz_1(i) = AE(vzz_cal, vzz_ref(i));
    ae_vxx_1(i) = AE(vxx_cal, vxx_ref(i));
    ae_vyy_1(i) = AE(vyy_cal, vyy_ref(i));
    ae_vzzz_1(i) = AE(vzzz_cal, vzzz_ref(i));
    ae_vxxz_1(i) = AE(vxxz_cal, vxxz_ref(i));
    ae_vyyz_1(i) = AE(vyyz_cal, vyyz_ref(i));
    L1_1(i) = log10(abs(vxx_cal + vyy_cal + vzz_cal));
    L2_1(i) = log10(abs(vzzz_cal + vxxz_cal + vyyz_cal));

    re_v_1(i) = RE(v_cal, v_ref(i));
    re_vz_1(i) = RE(vz_cal, vz_ref(i));
    re_vzz_1(i) = RE(vzz_cal, vzz_ref(i));
    re_vxx_1(i) = RE(vxx_cal, vxx_ref(i));
    re_vyy_1(i) = RE(vyy_cal, vyy_ref(i));
    re_vzzz_1(i) = RE(vzzz_cal, vzzz_ref(i));
    re_vxxz_1(i) = RE(vxxz_cal, vxxz_ref(i));
    re_vyyz_1(i) = RE(vyyz_cal, vyyz_ref(i));
end

% 1/12
dlat = 1 / 12;
nmax = 2160;
vnm = sha(lon0, lat0, nmax, dlat, dlat);

ae_v_12 = zeros(num, 1);
ae_vz_12 = zeros(num, 1);
ae_vzz_12 = zeros(num, 1);
ae_vxx_12 = zeros(num, 1);
ae_vyy_12 = zeros(num, 1);
ae_vzzz_12 = zeros(num, 1);
ae_vxxz_12 = zeros(num, 1);
ae_vyyz_12 = zeros(num, 1);
L1_12 = zeros(num, 1);
L2_12 = zeros(num, 1);

re_v_12 = zeros(num, 1);
re_vz_12 = zeros(num, 1);
re_vzz_12 = zeros(num, 1);
re_vxx_12 = zeros(num, 1);
re_vyy_12 = zeros(num, 1);
re_vzzz_12 = zeros(num, 1);
re_vxxz_12 = zeros(num, 1);
re_vyyz_12 = zeros(num, 1);

load ../data/vref_h_1_12.mat

for i = 1 : num
    h = H(i);
    v_cal = shs(vnm, nmax, h, 'v');
    vz_cal = shs(vnm, nmax, h, 'vz');
    vzz_cal = shs(vnm, nmax, h, 'vzz');
    vxx_cal = shs(vnm, nmax, h, 'vxx');
    vyy_cal = shs(vnm, nmax, h, 'vyy');
    vzzz_cal = shs(vnm, nmax, h, 'vzzz');
    vxxz_cal = shs(vnm, nmax, h, 'vxxz');
    vyyz_cal = shs(vnm, nmax, h, 'vyyz');

    ae_v_12(i) = AE(v_cal, v_ref(i));
    ae_vz_12(i) = AE(vz_cal, vz_ref(i));
    ae_vzz_12(i) = AE(vzz_cal, vzz_ref(i));
    ae_vxx_12(i) = AE(vxx_cal, vxx_ref(i));
    ae_vyy_12(i) = AE(vyy_cal, vyy_ref(i));
    ae_vzzz_12(i) = AE(vzzz_cal, vzzz_ref(i));
    ae_vxxz_12(i) = AE(vxxz_cal, vxxz_ref(i));
    ae_vyyz_12(i) = AE(vyyz_cal, vyyz_ref(i));
    L1_12(i) = log10(abs(vxx_cal + vyy_cal + vzz_cal));
    L2_12(i) = log10(abs(vzzz_cal + vxxz_cal + vyyz_cal));

    re_v_12(i) = RE(v_cal, v_ref(i));
    re_vz_12(i) = RE(vz_cal, vz_ref(i));
    re_vzz_12(i) = RE(vzz_cal, vzz_ref(i));
    re_vxx_12(i) = RE(vxx_cal, vxx_ref(i));
    re_vyy_12(i) = RE(vyy_cal, vyy_ref(i));
    re_vzzz_12(i) = RE(vzzz_cal, vzzz_ref(i));
    re_vxxz_12(i) = RE(vxxz_cal, vxxz_ref(i));
    re_vyyz_12(i) = RE(vyyz_cal, vyyz_ref(i));
end
H = H / 1e3;
%% plot
fontsize = 9;
fig = figure;
mksz = 18;
linewidth = 1;
outer = tiledlayout(fig, 30, 1, "Padding", "tight", 'TileSpacing', 'tight');
inner = tiledlayout(outer, 2, 2, "Padding", "tight", 'TileSpacing', 'tight');
inner.Layout.TileSpan = [29, 1];

nexttile(inner)
scatter(H, ae_v_1, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(H, ae_vz_1, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(H, ae_vzz_1, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(H, ae_vxx_1, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(H, ae_vyy_1, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(H, ae_vzzz_1, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(H, ae_vxxz_1, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(H, ae_vyyz_1, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
scatter(H, L1_1, mksz, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
scatter(H, L2_1, mksz, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
ylim([-40, 0]);
xlim([0, 500])
grid on
box on
title('(a)');
ylabel('Absolute errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
scatter(H, re_v_1, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(H, re_vz_1, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(H, re_vzz_1, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(H, re_vxx_1, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(H, re_vyy_1, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(H, re_vzzz_1, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(H, re_vxxz_1, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(H, re_vyyz_1, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
ylim([-15, 5]);
xlim([0, 500])

ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
grid on
box on
title('(b)');
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
scatter(H, ae_v_12, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(H, ae_vz_12, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(H, ae_vzz_12, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(H, ae_vxx_12, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(H, ae_vyy_12, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(H, ae_vzzz_12, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(H, ae_vxxz_12, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(H, ae_vyyz_12, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
scatter(H, L1_12, mksz, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
scatter(H, L2_12, mksz, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
ylim([-40, 0]);
xlim([0, 500])
grid on
box on
title('(c)');
ylabel('Absolute errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
scatter(H, re_v_12, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(H, re_vz_12, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(H, re_vzz_12, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(H, re_vxx_12, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(H, re_vyy_12, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(H, re_vzzz_12, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(H, re_vxxz_12, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(H, re_vyyz_12, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
ylim([-15, 5]);
xlim([0, 500])

ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
grid on
box on
title('(d)');
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);
xlabel(inner, 'Observation height $h$ (km)', 'Interpreter', 'latex', 'fontsize', fontsize);

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
scatter(nan, nan, mksz, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
axis off;
leg = legend('$V$', '$V_z$', '$V_{zz}$', '$V_{xx}$', '$V_{yy}$', '$V_{zzz}$', '$V_{xxz}$', '$V_{yyz}$', '$\delta L_1$','$\delta L_2$', 'fontsize', fontsize, 'NumColumns', 6, 'Interpreter', 'latex');
leg.Box = "off";
leg.Layout.Tile = 'South';
set(gcf,'Units','centimeters','Position', [1.3 3 13, 15]);
exportgraphics(gcf, '../manuscript/marked/tess_h.pdf');
%% nmax
h = 260 * 1e3;
dlat = 1;
lat1 = lat0 - dlat / 2; lat2 = lat0 + dlat / 2;

load ../data/vref_1_h260.mat

Nmax = 0 : 50 : 2000;
vnm = sha(lon0, lat0, Nmax(end), dlat, dlat);

num = length(Nmax);
ae_nmax_1_260_v = zeros(num, 1);
ae_nmax_1_260_vz = zeros(num, 1);
ae_nmax_1_260_vzz = zeros(num, 1);
ae_nmax_1_260_vxx = zeros(num, 1);
ae_nmax_1_260_vyy = zeros(num, 1);
ae_nmax_1_260_vzzz = zeros(num, 1);
ae_nmax_1_260_vxxz = zeros(num, 1);
ae_nmax_1_260_vyyz = zeros(num, 1);
L1_nmax_1_260 = zeros(num, 1);
L2_nmax_1_260 = zeros(num, 1);

re_nmax_1_260_v = zeros(num, 1);
re_nmax_1_260_vz = zeros(num, 1);
re_nmax_1_260_vzz = zeros(num, 1);
re_nmax_1_260_vxx = zeros(num, 1);
re_nmax_1_260_vyy = zeros(num, 1);
re_nmax_1_260_vzzz = zeros(num, 1);
re_nmax_1_260_vxxz = zeros(num, 1);
re_nmax_1_260_vyyz = zeros(num, 1);

for i = 1 : num
    nmax = Nmax(i);
    v_cal = shs(vnm, nmax, h, 'v');
    vz_cal = shs(vnm, nmax, h, 'vz');
    vzz_cal = shs(vnm, nmax, h, 'vzz');
    vxx_cal = shs(vnm, nmax, h, 'vxx');
    vyy_cal = shs(vnm, nmax, h, 'vyy');
    vzzz_cal = shs(vnm, nmax, h, 'vzzz');
    vxxz_cal = shs(vnm, nmax, h, 'vxxz');
    vyyz_cal = shs(vnm, nmax, h, 'vyyz');

    ae_nmax_1_260_v(i) = AE(v_cal, v_ref);
    ae_nmax_1_260_vz(i) = AE(vz_cal, vz_ref);
    ae_nmax_1_260_vzz(i) = AE(vzz_cal, vzz_ref);
    ae_nmax_1_260_vxx(i) = AE(vxx_cal, vxx_ref);
    ae_nmax_1_260_vyy(i) = AE(vyy_cal, vyy_ref);
    ae_nmax_1_260_vzzz(i) = AE(vzzz_cal, vzzz_ref);
    ae_nmax_1_260_vxxz(i) = AE(vxxz_cal, vxxz_ref);
    ae_nmax_1_260_vyyz(i) = AE(vyyz_cal, vyyz_ref);
    L1_nmax_1_260(i) = log10(abs(vxx_cal + vyy_cal + vzz_cal));
    L2_nmax_1_260(i) = log10(abs(vzzz_cal + vxxz_cal + vyyz_cal));

    re_nmax_1_260_v(i) = RE(v_cal, v_ref);
    re_nmax_1_260_vz(i) = RE(vz_cal, vz_ref);
    re_nmax_1_260_vzz(i) = RE(vzz_cal, vzz_ref);
    re_nmax_1_260_vxx(i) = RE(vxx_cal, vxx_ref);
    re_nmax_1_260_vyy(i) = RE(vyy_cal, vyy_ref);
    re_nmax_1_260_vzzz(i) = RE(vzzz_cal, vzzz_ref);
    re_nmax_1_260_vxxz(i) = RE(vxxz_cal, vxxz_ref);
    re_nmax_1_260_vyyz(i) = RE(vyyz_cal, vyyz_ref);
end

h = 10 * 1e3;
load ../data/vref_1_12_h10.mat

Nmax = 0 : 50 : 10000;
num = length(Nmax);
ae_nmax_12_10_v = zeros(num, 1);
ae_nmax_12_10_vz = zeros(num, 1);
ae_nmax_12_10_vzz = zeros(num, 1);
ae_nmax_12_10_vxx = zeros(num, 1);
ae_nmax_12_10_vyy = zeros(num, 1);
ae_nmax_12_10_vzzz = zeros(num, 1);
ae_nmax_12_10_vxxz = zeros(num, 1);
ae_nmax_12_10_vyyz = zeros(num, 1);
L1_nmax_12_10 = zeros(num, 1);
L2_nmax_12_10 = zeros(num, 1);

re_nmax_12_10_v = zeros(num, 1);
re_nmax_12_10_vz = zeros(num, 1);
re_nmax_12_10_vzz = zeros(num, 1);
re_nmax_12_10_vxx = zeros(num, 1);
re_nmax_12_10_vyy = zeros(num, 1);
re_nmax_12_10_vzzz = zeros(num, 1);
re_nmax_12_10_vxxz = zeros(num, 1);
re_nmax_12_10_vyyz = zeros(num, 1);

vnm = vnm_1_2;
for i = 1 : num
    nmax = Nmax(i);
    v_cal = shs(vnm, nmax, h, 'v');
    vz_cal = shs(vnm, nmax, h, 'vz');
    vzz_cal = shs(vnm, nmax, h, 'vzz');
    vxx_cal = shs(vnm, nmax, h, 'vxx');
    vyy_cal = shs(vnm, nmax, h, 'vyy');
    vzzz_cal = shs(vnm, nmax, h, 'vzzz');
    vxxz_cal = shs(vnm, nmax, h, 'vxxz');
    vyyz_cal = shs(vnm, nmax, h, 'vyyz');

    ae_nmax_12_10_v(i) = AE(v_cal, v_ref);
    ae_nmax_12_10_vz(i) = AE(vz_cal, vz_ref);
    ae_nmax_12_10_vzz(i) = AE(vzz_cal, vzz_ref);
    ae_nmax_12_10_vxx(i) = AE(vxx_cal, vxx_ref);
    ae_nmax_12_10_vyy(i) = AE(vyy_cal, vyy_ref);
    ae_nmax_12_10_vzzz(i) = AE(vzzz_cal, vzzz_ref);
    ae_nmax_12_10_vxxz(i) = AE(vxxz_cal, vxxz_ref);
    ae_nmax_12_10_vyyz(i) = AE(vyyz_cal, vyyz_ref);
    L1_nmax_12_10(i) = log10(abs(vxx_cal + vyy_cal + vzz_cal));
    L2_nmax_12_10(i) = log10(abs(vzzz_cal + vxxz_cal + vyyz_cal));

    re_nmax_12_10_v(i) = RE(v_cal, v_ref);
    re_nmax_12_10_vz(i) = RE(vz_cal, vz_ref);
    re_nmax_12_10_vzz(i) = RE(vzz_cal, vzz_ref);
    re_nmax_12_10_vxx(i) = RE(vxx_cal, vxx_ref);
    re_nmax_12_10_vyy(i) = RE(vyy_cal, vyy_ref);
    re_nmax_12_10_vzzz(i) = RE(vzzz_cal, vzzz_ref);
    re_nmax_12_10_vxxz(i) = RE(vxxz_cal, vxxz_ref);
    re_nmax_12_10_vyyz(i) = RE(vyyz_cal, vyyz_ref);
end
%% plot
fontsize = 9;
fig = figure;

outer = tiledlayout(fig, 30, 1, "Padding", "tight", 'TileSpacing', 'tight');
inner = tiledlayout(outer, 2, 2, "Padding", "tight", 'TileSpacing', 'tight');
inner.Layout.TileSpan = [29, 1];
mksz = 16;
linewidth = 0.8;
% a
Nmax = 0 : 50 : 2000;

nexttile(inner)
scatter(Nmax, ae_nmax_1_260_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(Nmax, ae_nmax_1_260_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_1_260_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_1_260_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_1_260_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_1_260_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_1_260_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_1_260_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
scatter(Nmax, L1_nmax_1_260, mksz, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
scatter(Nmax, L2_nmax_1_260, mksz, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);

xlabel(inner, 'Truncated harmonic degree $n_{\mathrm{max}}$', 'Interpreter', 'latex', 'fontsize', fontsize);
xlim([0, 2000]);
ylim([-40, 10]);
grid on
box on
title('(a)');
ylabel('Absolute errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

% b
nexttile(inner)
scatter(Nmax, re_nmax_1_260_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(Nmax, re_nmax_1_260_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_1_260_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_1_260_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_1_260_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_1_260_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_1_260_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_1_260_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);

ylim([-15, 10]);
xlim([0, 2000]);
ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
grid on
box on
title('(b)');
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);


mksz = 8;
linewidth = 0.5;
Nmax = 0 : 50 : 10000;
% c
nexttile(inner)
scatter(Nmax, ae_nmax_12_10_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(Nmax, ae_nmax_12_10_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_12_10_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_12_10_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_12_10_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_12_10_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_12_10_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(Nmax, ae_nmax_12_10_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
scatter(Nmax, L1_nmax_12_10, mksz, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
scatter(Nmax, L2_nmax_12_10, mksz, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
xlim([0, 10000]);
xticks(0 : 2500 : 10000);
ylim([-40, 10]);
grid on
box on
title('(c)');
ylabel('Absolute errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

% d
nexttile(inner)
scatter(Nmax, re_nmax_12_10_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(Nmax, re_nmax_12_10_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_12_10_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_12_10_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_12_10_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_12_10_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_12_10_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(Nmax, re_nmax_12_10_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);

ylim([-15, 10]);
xlim([0, 10000]);
xticks(0 : 2500 : 10000);
ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
grid on
box on
title('(d)');
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
scatter(nan, nan, mksz, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);
scatter(nan, nan, mksz, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', linewidth);

axis off;
leg = legend('$V$', '$V_z$', '$V_{zz}$', '$V_{xx}$', '$V_{yy}$', '$V_{zzz}$', '$V_{xxz}$', '$V_{yyz}$', '$\delta L_1$','$\delta L_2$', 'fontsize', fontsize, 'NumColumns', 6, 'Interpreter', 'latex');
leg.Box = "off";
leg.Layout.Tile = 'South';

set(gcf,'Units','centimeters','Position', [1.3, 3, 14, 15]);
exportgraphics(gcf, '../manuscript/marked/tess_nmax.pdf');
%% theta0
nmax = 1000;
h = 260 * 1e3;
dlat = 1 / 12;
Lat0 = 0 : dlat : 90;
Lat0(end) = 90 - dlat / 2;
lon0 = 0;
num = length(Lat0);
ae_lat0_v = zeros(num, 1);
ae_lat0_vz = zeros(num, 1);
ae_lat0_vzz = zeros(num, 1);
ae_lat0_vxx = zeros(num, 1);
ae_lat0_vyy = zeros(num, 1);
ae_lat0_vzzz = zeros(num, 1);
ae_lat0_vxxz = zeros(num, 1);
ae_lat0_vyyz = zeros(num, 1);
L1_lat0= zeros(num, 1);
L2_lat0 = zeros(num, 1);

re_lat0_v = zeros(num, 1);
re_lat0_vz = zeros(num, 1);
re_lat0_vzz = zeros(num, 1);
re_lat0_vxx = zeros(num, 1);
re_lat0_vyy = zeros(num, 1);
re_lat0_vzzz = zeros(num, 1);
re_lat0_vxxz = zeros(num, 1);
re_lat0_vyyz = zeros(num, 1);

load ../data/vref_lat0.mat

for i = 1 : num
    lat0 = Lat0(i);
    vnm = sha(lon0, lat0, nmax, dlat, dlat);
    v_cal = shs(vnm, nmax, h, 'v');
    vz_cal = shs(vnm, nmax, h, 'vz');
    vzz_cal = shs(vnm, nmax, h, 'vzz');
    vxx_cal = shs(vnm, nmax, h, 'vxx');
    vyy_cal = shs(vnm, nmax, h, 'vyy');
    vzzz_cal = shs(vnm, nmax, h, 'vzzz');
    vxxz_cal = shs(vnm, nmax, h, 'vxxz');
    vyyz_cal = shs(vnm, nmax, h, 'vyyz');

    ae_lat0_v(i) = AE(v_cal, v_ref(i));
    ae_lat0_vz(i) = AE(vz_cal, vz_ref(i));
    ae_lat0_vzz(i) = AE(vzz_cal, vzz_ref(i));
    ae_lat0_vxx(i) = AE(vxx_cal, vxx_ref(i));
    ae_lat0_vyy(i) = AE(vyy_cal, vyy_ref(i));
    ae_lat0_vzzz(i) = AE(vzzz_cal, vzzz_ref(i));
    ae_lat0_vxxz(i) = AE(vxxz_cal, vxxz_ref(i));
    ae_lat0_vyyz(i) = AE(vyyz_cal, vyyz_ref(i));
    L1_lat0(i) = log10(abs(vxx_cal + vyy_cal + vzz_cal));
    L2_lat0(i) = log10(abs(vzzz_cal + vxxz_cal + vyyz_cal));

    re_lat0_v(i) = RE(v_cal, v_ref(i));
    re_lat0_vz(i) = RE(vz_cal, vz_ref(i));
    re_lat0_vzz(i) = RE(vzz_cal, vzz_ref(i));
    re_lat0_vxx(i) = RE(vxx_cal, vxx_ref(i));
    re_lat0_vyy(i) = RE(vyy_cal, vyy_ref(i));
    re_lat0_vzzz(i) = RE(vzzz_cal, vzzz_ref(i));
    re_lat0_vxxz(i) = RE(vxxz_cal, vxxz_ref(i));
    re_lat0_vyyz(i) = RE(vyyz_cal, vyyz_ref(i));
end
CLat0 = 90 - Lat0;
%%
fontsize = 9;
fig = figure;
outer = tiledlayout(fig, 30, 1, "Padding", "tight", 'TileSpacing', 'tight');
inner = tiledlayout(outer, 3, 1, "Padding", "tight", 'TileSpacing', 'tight');
inner.Layout.TileSpan = [29, 1];
mksz = 4;
linewidth = 0.5;
% a
nexttile(inner)
plot(CLat0, log10(abs(v_ref)), 'Color', cmp(1, :), 'LineWidth', 1.2);
hold on
plot(CLat0, log10(abs(vz_ref)), 'Color', cmp(2, :), 'LineWidth', 1.2);
plot(CLat0, log10(abs(vzz_ref)), 'Color', cmp(3, :), 'LineWidth', 1.2);
plot(CLat0, log10(abs(vxx_ref)), 'Color',  cmp(4, :), 'LineWidth', 1.2);
plot(CLat0, log10(abs(vyy_ref)), 'Color', cmp(5, :), 'LineWidth', 1.2);
plot(CLat0, log10(abs(vzzz_ref)), 'Color', cmp(6, :), 'LineWidth', 1.2);
plot(CLat0, log10(abs(vxxz_ref)), 'Color', cmp(7, :), 'LineWidth', 1.2);
plot(CLat0, log10(abs(vyyz_ref)), 'Color', cmp(8, :), 'LineWidth', 1.2);
leg = legend('$V$', '$V_z$', '$V_{zz}$', '$V_{xx}$', '$V_{yy}$', '$V_{zzz}$', '$V_{xxz}$', '$V_{yyz}$', 'fontsize', fontsize, 'NumColumns', 10, 'Interpreter', 'latex', 'location', 'southeast');
leg.ItemTokenSize = [10, 10];

xlim([0, 90]);
ylim([-30, 0]);
xticks(0:15:90);
set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabel', {'$0^\circ$', '$15^\circ$', '$30^\circ$', '$45^\circ$', '$60^\circ$', '$75^\circ$', '$90^\circ$'})
ylabel('Absolute field values in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
grid on
box on
title('(a)');
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile(inner)
scatter(CLat0, ae_lat0_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(CLat0, ae_lat0_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(CLat0, ae_lat0_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(CLat0, ae_lat0_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(CLat0, ae_lat0_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(CLat0, ae_lat0_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(CLat0, ae_lat0_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(CLat0, ae_lat0_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);
scatter(CLat0, L1_lat0, mksz, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', 0.5);
scatter(CLat0, L2_lat0, mksz, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', 0.5);
xlabel(inner, 'Tesseroid central colatitude $\theta_0^\prime$', 'Interpreter', 'latex', 'fontsize', fontsize)
xlim([0, 90]);
xticks(0:15:90);
set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabel', {'$0^\circ$', '$15^\circ$', '$30^\circ$', '$45^\circ$', '$60^\circ$', '$75^\circ$', '$90^\circ$'})
ylim([-40, -10]);
grid on
box on
title('(b)');
ylabel('Absolute errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

% b
nexttile(inner)
scatter(CLat0, re_lat0_v, mksz, 'filled', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(1, :), 'LineWidth', linewidth);
hold on
scatter(CLat0, re_lat0_vz, mksz, 'o', 'MarkerEdgeColor', cmp(2, :), 'LineWidth', linewidth);
scatter(CLat0, re_lat0_vzz, mksz, 'filled', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(3, :), 'LineWidth', linewidth);
scatter(CLat0, re_lat0_vxx, mksz, 's', 'MarkerEdgeColor', cmp(4, :), 'LineWidth', linewidth);
scatter(CLat0, re_lat0_vyy, mksz, 'filled', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(5, :), 'LineWidth', linewidth);
scatter(CLat0, re_lat0_vzzz, mksz, '^', 'MarkerEdgeColor', cmp(6, :), 'LineWidth', linewidth);
scatter(CLat0, re_lat0_vxxz, mksz, 'filled', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmp(7, :), 'LineWidth', linewidth);
scatter(CLat0, re_lat0_vyyz, mksz, 'v', 'MarkerEdgeColor', cmp(8, :), 'LineWidth', linewidth);

xlim([0, 90]);
xticks(0:15:90);
ylim([-16,-4])
set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabel', {'$0^\circ$', '$15^\circ$', '$30^\circ$', '$45^\circ$', '$60^\circ$', '$75^\circ$', '$90^\circ$'})
ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)
grid on
box on
title('(c)');
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
scatter(nan, nan, mksz, 'x', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', 0.5);
scatter(nan, nan, mksz, '+', 'MarkerEdgeColor', cmp(9, :), 'LineWidth', 0.5);

axis off;
leg = legend('$V$', '$V_z$', '$V_{zz}$', '$V_{xx}$', '$V_{yy}$', '$V_{zzz}$', '$V_{xxz}$', '$V_{yyz}$', '$\delta L_1$','$\delta L_2$', 'fontsize', fontsize, 'NumColumns', 6, 'Interpreter', 'latex');
leg.Box = "off";
leg.Layout.Tile = 'South';

set(gcf,'Units','centimeters','Position', [1.3 3 13, 17]);
exportgraphics(gcf, '../manuscript/marked/tess_lat0.pdf');