function rs = power_log10(SHCs, GM, r, a, type)
%{
    Used for calculating the normalized signal power spectra
Input:
    SHCs: [n, m, Cnm, Snm]
    GM: G*M constants
    r: the distance from the center at the point [m]
    a: semi-major axis of the ellipsoid [m]
    type: calculation type
Return:
    rs: (:, 1) the normalized signal power spectra, sort by degree
%}

arguments
    SHCs (:, :)
    GM (1, 1)
    r (1, 1)
    a (1, 1)
    type
end

nmax = SHCs(end, 1);
SHCs = cnm2sc(SHCs, nmax);

n = (0 : nmax)';
sn = sum(SHCs .^ 2, 2);
digits(100);
a = vpa(a);
r = vpa(r);

switch type
    case 'v'
        rs = sn .* (a / r) .^ (2 * n + 2) * (GM / a) ^ 2;
    case 'vz'
        rs = sn .* (a / r) .^ (2 * n + 4) .* (n + 1) .^ 2 * (GM / a ^ 2) ^ 2;
    case 'vzz'
        rs = sn .*  (a / r) .^ (2 * n + 4) .* (n + 2) .^ 2 .* (n + 1) .^ 2 * (GM / a ^ 3) ^ 2;
    case 'vzzz'
        rs = sn .* (a / r) .^ (2 * n + 8) .* (n + 3) .^ 2 .* (n + 2) .^ 2 .* (n + 1) .^ 2 * (GM / a ^ 4) ^ 2;
end
rs = rs / sum(rs);
rs = double(log10(rs));

end

