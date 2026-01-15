clear
close all
clc
%% parameter
a = 6371e3;
M = 5.965e24;
G = 6.672e-11;
h = 10e3;
dlon = 2;
dlat = dlon;
lon = -180 + dlon / 2 : dlon : 180;
lat = 90 - dlat / 2: -dlat : -90;
[lon_gd, lat_gd] = meshgrid(lon, lat);
nmax = 20000;
%% data
load('../data/ETOPO2022_topo_blockmean_to_2d.mat');

r2 = a + topo;
r1 = r2 * 0 + a;
rho = r2 * 0 + 2670;
ind = topo < 0;
rho(ind) = -1640;

clat = linspace(-89, 89, 8); clon = linspace(-179, 179, 16);
[clon_gd, clat_gd] = meshgrid(clon, clat);

Phis = lat_gd(:) - dlat / 2;
Phin = lat_gd(:) + dlat / 2;
Lamw = lon_gd(:) - dlon / 2;
Lame = lon_gd(:) + dlon / 2;

HT = r2 - a;
HB = r1 - a;
Rho = rho;

t = HT;
HT(ind) = HB(ind);
HB(ind) = t(ind);
Rho(ind) = -Rho(ind);

Mod = [Phis, Phin, Lamw, Lame, HB(:), HT(:), Rho(:)];
Obs = [clat_gd(:), clon_gd(:), clon_gd(:) * 0 + h];
nmods = size(Mod, 1);
nobs = size(Obs, 1);
%% write input file
fname = 'input_xtessgc.dat';
fid = fopen(fname, 'w');
assert(fid > 0, 'can not open：%s', fname);

fprintf(fid, '%.15E %.15E\n', G, a);
fprintf(fid, '%d %d\n\n', nmods, nobs);

for i = 1:nmods
    fprintf(fid, '%.15E %.15E %.15E %.15E %.15E %.15E %.15E\n', Mod(i, :));
end
fprintf(fid, '\n\n');
% computation points
for i = 1:nobs
    fprintf(fid, '%.15E %.15E %.15E\n', Obs(i, :));
end
fclose(fid);
%% read data
spa_field = load('../data/field_derived_by_xtessgcf90.dat');
spa_field(:, [3, 6, 9, 12, 14, 16, 19]) = -spa_field(:, [3, 6, 9, 12, 14, 16, 19]); % 统一方向
% V, Vx, Vy, Vz, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz, Vxxx, Vxxy, Vxxz, Vxyz, Vxyy, Vyyy, Vyyz, Vxzz, Vyzz, Vzzz
REs = @(x1, x2) log10(abs( (x1 - x2) ./ x2) ) ;
%%
spe_field = zeros(size(spa_field));

tic
vnm = tess2Vnm(lon, lat, r1, r2, rho, nmax, "a", a, "M", M);
toc
shs = @(xtype) reshape(shs_grid(vnm, clon, clat, nmax, 0, "height", h, "fldTyp", xtype, "GM", G * M, "a", a), [], 1);

spe_field(:, 1) = shs('v');
spe_field(:, 2) = shs('vx');
spe_field(:, 3) = shs('vy');
spe_field(:, 4) = shs('vz');
spe_field(:, 5) = shs('vxx');
spe_field(:, 6) = shs('vxy');
spe_field(:, 7) = shs('vxz');
spe_field(:, 8) = shs('vyy');
spe_field(:, 9) = shs('vyz');
spe_field(:, 10) = shs('vzz');
spe_field(:, 11) = shs('vxxx');
spe_field(:, 12) = shs('vxxy');
spe_field(:, 13) = shs('vxxz');
spe_field(:, 14) = shs('vxyz');
spe_field(:, 15) = shs('vxyy');
spe_field(:, 16) = shs('vyyy');
spe_field(:, 17) = shs('vyyz');
spe_field(:, 18) = shs('vxzz');
spe_field(:, 19) = shs('vyzz');
spe_field(:, 20) = shs('vzzz');

re_field = REs(spe_field, spa_field);


meanVals = nan(1, 20);
minVals = meanVals;
maxVals = meanVals;

for i = 1 : 20
    t = re_field(:, i);
    t(isinf(t)) = [];
    meanVals(i) = mean(t);
    minVals(i) = min(t);
    maxVals(i) = max(t);
end
%%
fontsize = 9;
nComp = 20;

labels = {'$V$', ...
    '$V_x$','$V_y$','$V_z$', ...
    '$V_{xx}$','$V_{xy}$','$V_{xz}$','$V_{yy}$','$V_{yz}$','$V_{zz}$', ...
    '$V_{xxx}$','$V_{xxy}$','$V_{xxz}$','$V_{xyy}$','$V_{xyz}$', ...
    '$V_{xzz}$','$V_{yyy}$','$V_{yyz}$','$V_{yzz}$','$V_{zzz}$'};

fig = figure;
tiledlayout(1, 1, "Padding", "compact", 'TileSpacing', 'compact', Parent=fig);
nexttile
x = 1 : nComp;

yl = [-16, 0];

ms = 4; lw = 1;
hold on
c1 = [0.87 0.92 1.00];  % blue-ish
c2 = [0.90 1.00 0.90];  % green-ish
c3 = [1.00 0.90 0.90];  % red-ish
patch([0.5  4.5  4.5  0.5], [yl(1) yl(1) yl(2) yl(2)], c1, 'EdgeColor','none','FaceAlpha',0.50);
patch([4.5 10.5 10.5 4.5 ], [yl(1) yl(1) yl(2) yl(2)], c2, 'EdgeColor','none','FaceAlpha',0.50);
patch([10.5 20.5 20.5 10.5], [yl(1) yl(1) yl(2) yl(2)], c3, 'EdgeColor','none','FaceAlpha',0.50);

errorbar(x, meanVals, abs(meanVals - minVals), abs(maxVals - meanVals), 'o', "MarkerSize", ms,...
    "MarkerEdgeColor","none", "MarkerFaceColor", [0.2, 0.2, 0.2], 'LineWidth', 0.9, 'Color', [0.2, 0.2, 0.2])
xlim([0, 21]);
ylim([-16, 1]);
yticks(-16 : 4 : 0);
set(gca, 'XTick', x, 'XTickLabel', labels, ...
    'TickLabelInterpreter','latex', 'FontSize', fontsize, 'GridLineStyle',':', 'GridAlpha', 0.5, ...
    'Box','off');
set(gca,'XTickLabelRotation', 90);
grid on
ylabel('Relative errors in $\log_{10}$', 'Interpreter', 'latex', 'fontsize', fontsize)

txtY = yl(2) - 1.1;

set(gcf, 'Units', 'centimeters', 'Position', [1.3, 3, 14, 8]);

exportgraphics(gcf, '../manuscript/marked/relief_REs.pdf');

