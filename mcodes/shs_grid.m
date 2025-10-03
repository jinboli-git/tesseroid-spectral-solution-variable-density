function fld = shs_grid(SHCs, lon, lat, nmax, nmin, parm)
%{
    Calculate the spherical harmonic synthesis of any rectangular grid
Input:
   SHCs: Spherical harmonic coefficient, [n, m, Cnm, Snm]
   lon: the longitude vector of the computation grid, [deg]
   lat: the latitude vector of the computation grid, [deg]
   nmax: maximum computational degree/order
   nmin: minimum computational degree/order
   parm
    ::
   fldTyp:  the type of calculation
   'shs', [none] only spherical harmonic synthesis
   'v', [m^2 s^-2] potiental
   'vi', first-order gradients [m s^-2]
   'vij', second-order gradients [s^-2]
   'vijk', third-order gradients [m^-1 s^-2]
    {i, j, k} = {x, y, z}
    a: radius of the reference sphere [m]
    GM: G*M constant
    height: computation height [m]
Ref:
    Sneeuw N (1994) Global spherical harmonic analysis by least-squares and numerical quadrature methods in historical perspective. Geophys J Int 118:707–716. https://doi.org/10.1111/j.1365-246X.1994.tb03995.x
    Eshagh M (2008) Non-singular expressions for the vector and the gradient tensor of gravitation in a geocentric spherical frame. Comput Geosci 34:1762–1768. https://doi.org/10.1016/j.cageo.2008.02.022
    Petrovskaya MS, Vershkov AN (2006) Non-Singular Expressions for the Gravity Gradients in the Local North-Oriented and Orbital Reference Frames. J Geodesy 80:117–127. https://doi.org/10.1007/s00190-006-0031-2
    Hamáčková E, Šprlák M, Pitoňák M, Novák P (2016) Non-singular expressions for the spherical harmonic synthesis of gravitational curvatures in a local north-oriented reference frame. Comput Geosci 88:152–162. https://doi.org/10.1016/j.cageo.2015.12.011
%}

arguments
    SHCs (:, :)
    lon (1, :)
    lat (:, 1)
    nmax = []
    nmin = 0;
    parm.fldTyp char {mustBeMember(parm.fldTyp, {'shs', 'v', 'vx', 'vy', 'vz', 'vxx', 'vxy', 'vxz', 'vyy', 'vyz', 'vzz', 'vxxx', 'vxxy', 'vxxz', 'vyyy', 'vxyy', 'vyyz', 'vxxz', 'vxyz', 'vyyz', 'vyzz', 'vxzz', 'vzzz'})} = 'shs'
    parm.a (1, 1) = 6371 * 1e3
    parm.GM (1, 1) = 5.965e24 * 6.672e-11
    parm.height (1, 1) = 0
end

nthe = length(lat);
r = parm.a + parm.height;

n = SHCs(:, 1);
nmaxSHCs = n(end);

if isempty(nmax)
    nmax = nmaxSHCs;
end
if nmin > nmax
    nmin = nmax;
end
if nmax > nmaxSHCs
    nmax = nmaxSHCs;
end

SHCs = SHCs(n >= nmin & n <= nmax, :);
n = SHCs(:, 1);
m = SHCs(:, 2);
listind = (1 : length(n))';
SHCs = SHCs(:, 3 : 4);

costhe = cosd(90 - lat);
Pnm = @(x) Pnm_Bel(nmax, costhe(x));
u = kapa(n, parm.fldTyp, parm.GM, parm.a, r);

SHCs = SHCs .* u;
SHCs = complex(SHCs(:, 1), SHCs(:, 2));

AB = zeros(nthe, nmax + 1, 'like', 1 + 1i);
mlam = lon .* (0 : nmax)';

% 1st step
switch parm.fldTyp
    case {'shs', 'v', 'vz', 'vzz', 'vzzz'}
        for ii = 1 : nthe
            pnm = Pnm(ii) .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 0 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vx'
        i = m == 0; ni = n(i);
        anm = sqrt( (1 + (m == 1)) .* (n + m) .* (n - m + 1) ) / 2;
        anm(i) = 0;
        bnm = -sqrt( (n - m) .* (n + m + 1) ) / 2;
        bnm(i) = -sqrt(ni .* (ni + 1) / 2);
        SHCs = -SHCs;
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = anm .* circshift(pnm, 1) + bnm .* circshift(pnm, -1);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 0 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vy'
        i = m == 0;
        t = sqrt((2 * n + 1) ./ (2 * n - 1)) / 2;
        anm = t .* sqrt( (1 + (m == 1)) .* (n + m) .* (n + m - 1) );
        bnm = t .* sqrt( (n - m) .* (n - m - 1) );
        anm(i) = 0; bnm(i) = 0;
        id = listind - n;
        SHCs = SHCs * 1i;
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = pnm(id);
            pnm = anm .* circshift(pnm, 1) + bnm .* circshift(pnm, -1);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 1 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vxx'
        i = m < 2; j = m == 0;
        ni = n(i); mi = m(i);
        anm = sqrt( (1 + (m == 2)) .* (n .^ 2 - (m - 1) .^ 2) .* (n + m) .* (n - m + 2) ) / 4;
        anm(i) = 0;
        bnm = (n .^ 2 + m .^ 2 + 3 * n + 2) / 2;
        bnm(i) = (ni + mi + 1) .* (ni + mi + 2) ./ (mi + 1) / 2;
        bnm = bnm - (n + 1) .* (n + 2);
        cnm = sqrt( (n .^ 2 - (m + 1) .^ 2) .* (n - m) .* (n + m + 2) ) / 4;
        cnm(j) = cnm(j) * sqrt(2);
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = anm .* circshift(pnm, 2) + bnm .* pnm + cnm .* circshift(pnm, -2);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 0 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vxy'
        t = sqrt((2 * n + 1) ./ (2 * n - 1));
        i0 = m == 0; i = m == 1; ni = n(i); ti = t(i);
        anm = -t .* sqrt((1 + (m == 2)) .* (n .^ 2 - (m - 1) .^ 2) .* (n + m) .* (n + m - 2) ) / 4;
        anm(i) = 0;
        bnm = m .* t .* sqrt((n + m) .* (n - m)) / 2;
        bnm(i) = ti .* sqrt((ni + 1) .* (ni - 1)) .* (ni + 2) / 4;
        cnm = t .* sqrt( (n .^ 2 - (m + 1) .^ 2) .* (n - m) .* (n - m - 2)) / 4;
        cnm(i) = ti .* sqrt( (ni - 3) .* (ni - 2) .* (ni - 1) .* (ni + 2) ) / 4;
        anm(i0) = 0; bnm(i0) = 0; cnm(i0) = 0;
        id = listind - n;
        SHCs = SHCs * 1i;
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = pnm(id);
            pnm = anm .* circshift(pnm, 2) + bnm .* pnm + cnm .* circshift(pnm, -2);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 1 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vxz'
        i = m == 0; ni = n(i);
        anm = (n + 2) .* sqrt( (1 + (m == 1)) .* (n + m) .* (n - m + 1) ) / 2;
        anm(i) = 0;
        bnm = -(n + 2) .* sqrt( (n - m) .* (n + m + 1) ) / 2;
        bnm(i) = -(ni + 2) .* sqrt(ni .* (ni + 1) / 2);
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = anm .* circshift(pnm, 1) + bnm .* circshift(pnm, -1);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 0 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vyy'
        i = m < 2; j = m == 0;
        ni = n(i); mi = m(i);
        anm = sqrt( (1 + (m == 2)) .* (n .^ 2 - (m - 1) .^ 2) .* (n + m) .* (n - m + 2) ) / 4;
        anm(i) = 0;
        bnm = (n .^ 2 + m .^ 2 + 3 * n + 2) / 2;
        bnm(i) = (ni + mi + 1) .* (ni + mi + 2) ./ (mi + 1) / 2;
        cnm = sqrt( (n .^ 2 - (m + 1) .^ 2) .* (n - m) .* (n + m + 2) ) / 4;
        cnm(j) = cnm(j) * sqrt(2);
        SHCs = -SHCs;
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = anm .* circshift(pnm, 2) + bnm .* pnm + cnm .* circshift(pnm, -2);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 0 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vyz'
        i = m == 0;
        t = -(n + 2) .* sqrt((2 * n + 1) ./ (2 * n - 1)) / 2;
        anm = t .* sqrt( (1 + (m == 1)) .* (n + m) .* (n + m - 1) );
        bnm = t .* sqrt( (n - m) .* (n - m - 1) );
        anm(i) = 0; bnm(i) = 0;
        id = listind - n;
        SHCs = SHCs * 1i;
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = pnm(id);
            pnm = anm .* circshift(pnm, 1) + bnm .* circshift(pnm, -1);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 1 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vxxx'
        i = m < 3; i0 = m == 0; i1 = m == 1; i2 = m == 2;
        n0 = n(i0); n1 = n(i1); n2 = n(i2);
        anm = 1 / 8 * sqrt( (1 + (m == 3)) .* (n + m) .* (n + m - 1) .* (n + m - 2) .* (n - m + 1) .* (n - m + 2) .* (n - m + 3) );
        anm(i) = 0;
        bnm = 3 / 8 * sqrt( (n + m) .* (n - m + 1) ) .* (n + m + 2) .* (m - n - 3);
        bnm(i0) = 0;
        bnm(i1) = -3 / 8 * sqrt( 2 * n1 .* (n1 + 1) ) .* (n1 + 2) .* (n1 + 3);
        bnm(i2) = -1 / 2 * sqrt( (n2 - 1) .* (n2 + 2) ) .* (n2 + 1) .* (n2 + 3);
        cnm = -3 / 8 * sqrt( (n - m) .* (n + m + 1) ) .* (m - n - 2) .* (m + n + 3);
        cnm(i0) = 3 / 4 * sqrt(n0 .* (n0 + 1) / 2) .* (n0 + 2) .* (n0 + 3);
        cnm(i1) = 1 / 2 * sqrt( (n1 - 1) .* (n1 + 2) ) .* (n1 + 1) .* (n1 + 3);
        cnm(i2) = 3 / 8 * sqrt( (n2 - 2) .* (n2 + 3) ) .* n2 .* (n2 + 5);
        dnm = -1 / 8 * sqrt( (n - m) .* (n - m - 1) .* (n - m - 2) .* (n + m + 1) .* (n + m + 2) .* (n + m + 3) );
        dnm(i0) = -1 / 4 * sqrt( n0 .* (n0 - 1) .* (n0 - 2) .* (n0 + 1) .* (n0 + 2) .* (n0 + 3) / 2);
        SHCs = -SHCs;
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = anm .* circshift(pnm, 3) + bnm .* circshift(pnm, 1) + cnm .* circshift(pnm, -1) + dnm .* circshift(pnm, -3);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 0 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vxxy'
        i0 = m == 0; i1 = m == 1; i2 = m == 2;
        n1 = n(i1); n2 = n(i2);
        anm = -1 / 8 ./ m .* sqrt( (2 * n + 1) ./ (2 * n - 1) .* (1 + (m == 3)) .* (m + n) .* (n + m - 1) .* (n + m - 2) .* (n + m - 3) .* (n - m + 1) .* (n - m + 2) );
        anm(i0 | i1 | i2) = 0;
        bnm = -1 / 8 ./ m .* sqrt( (2 * n + 1) ./ (2 * n - 1) .* (m + n) .* (n + m - 1) ) .* (3 * m .* m - m .* (2 * n + 3) - n .* n -3 * (n + 2));
        bnm(i0) = 0;
        bnm(i1) = 1 / 8 * sqrt( 2 * n1 .* (n1 + 1) .* (2 * n1 + 1) ./ (2 * n1 - 1) ) .* (n1 + 2) .* (n1 + 3);
        bnm(i2) = 1 / 8 * sqrt( (n2 + 1) .* (n2 + 2) .* (2 * n2 + 1) ./ (2 * n2 - 1) ) .* n2 .* (n2 + 3);
        cnm = -1 / 8 ./ m .* sqrt( (2 * n + 1) ./ (2 * n - 1) .* (m - n) .* (m - n + 1) ) .* (3 * m .* m + m .* (2 * n + 3) - n .* n - 3 * (n + 2));
        cnm(i0 | i1) = 0;
        cnm(i2) = 1 / 16 * sqrt( (n2 - 2) .* (n2 - 3) .* (2 * n2 + 1) ./ (2 * n2 - 1) ) .* (n2 + 3) .* (n2 - 4);
        dnm = -1 / 8 ./ m .* sqrt( (2 * n + 1) ./ (2 * n - 1) .* (n - m) .* (n + m + 1) .* (n + m + 2) .* (n - m - 1) .* (n - m - 2) .* (n - m - 3) );
        dnm(i0) = 0;
        id = listind - n;
        SHCs = SHCs .* m * -1i;
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = pnm(id);
            pnm = anm .* circshift(pnm, 3) + bnm .* circshift(pnm, 1) + cnm .* circshift(pnm, -1) + dnm .* circshift(pnm, -3);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 1 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vxxz'
        i = m < 2; i0 = m == 0; i1 = m == 1;
        n0 = n(i0); n1 = n(i1);
        anm = -1 / 4 * sqrt( (1 + (m == 2)) .* (n + m) .* (n + m - 1) .* (n - m + 1) .* (n - m + 2) );
        anm(i) = 0;
        bnm = 1 / 2 * ((n + 1) .* (n + 2) - m .* m);
        bnm(i1) = 1 / 4 * (n1 + 2) .* (3 * n1 + 1);
        cnm = -1 / 4 * sqrt( (n - m) .* (n - m - 1) .* (n + m + 1) .* (n + m + 2) );
        cnm(i0) = -1 / 4 * sqrt( 2 * (n0 - 1) .* n0 .* (n0 + 1) .* (n0 + 2) );
        SHCs = SHCs .* (n + 3);
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = anm .* circshift(pnm, 2) + bnm .* pnm + cnm .* circshift(pnm, -2);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 0 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vxyy'
        i = m < 3; i0 = m == 0; i1 = m == 1; i2 = m == 2;
        n0 = n(i0); n1 = n(i1); n2 = n(i2);
        anm = 1 / 8 * sqrt( (1 + (m == 3)) .* (n + m) .* (n + m - 1) .* (n + m - 2) .* (n - m + 1) .* (n - m + 2) .* (n - m + 3) );
        anm(i) = 0;
        bnm = 1 / 8 * sqrt( (1 + (m == 1)) .* (n + m) .* (n - m + 1) ) .* ( 3 * m .* (m - 1) + (n + 2) .* (n + 3) );
        bnm(i0) = 0;
        bnm(i2) = 1 / 2 * sqrt( (n2 - 1) .* (n2 + 2) ) .* (n2 + 3);
        cnm = -1 / 8 * sqrt( (n - m) .* (n + m + 1) ) .* ( 3 * m .* (m + 1) + (n + 2) .* (n + 3) );
        cnm(i0) = -1 / 4 * sqrt(n0 .* (n0 + 1) / 2) .* (n0 + 2) .* (n0 + 3);
        cnm(i1) = -1 / 2 * sqrt( (n1 - 1) .* (n1 + 2) ) .* (n1 + 3);
        cnm(i2) = -1 / 8 * sqrt( (n2 - 2) .* (n2 + 3) ) .* (n2 .* n2 + 5 * n2 + 24);
        dnm = -1 / 8 * sqrt( (n - m) .* (n - m - 1) .* (n - m - 2) .* (n + m + 1) .* (n + m + 2) .* (n + m + 3) );
        dnm(i0) = -1 / 4 * sqrt( n0 .* (n0 - 1) .* (n0 - 2) .* (n0 + 1) .* (n0 + 2) .* (n0 + 3) / 2);
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = anm .* circshift(pnm, 3) + bnm .* circshift(pnm, 1) + cnm .* circshift(pnm, -1) + dnm .* circshift(pnm, -3);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 0 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vxyz'
        i = m < 2; i0 = m == 0; i1 = m == 1;
        n1 = n(i1);
        anm = -1 / 4 * sqrt( (2 * n + 1) ./ (2 * n - 1) .* (1 + (m == 2)) .* (n + m) .* (n + m - 1) .* (n + m - 2) .* (n - m + 1) );
        anm(i) = 0;
        bnm = m / 2 .* sqrt( (2 * n + 1) ./ (2 * n - 1) .* (n + m) .* (n - m) );
        bnm(i0) = 0;
        bnm(i1) = 1 / 4 * sqrt( (2 * n1 + 1) ./ (2 * n1 - 1) .* (n1 + 1) .* (n1 - 1) ) .* (n1 + 2);
        cnm = 1 / 4 * sqrt( (2 * n + 1) ./ (2 * n - 1) .* (n - m) .* (n - m - 1) .* (n - m - 2) .* (n + m + 1) );
        cnm(i0) = 0;
        id = listind - n;
        SHCs = SHCs .* (n + 3) * -1i;
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = pnm(id);
            pnm = anm .* circshift(pnm, 2) + bnm .* pnm + cnm .* circshift(pnm, -2);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 1 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vxzz'
        i = m == 0;
        n0 = n(i);
        anm = -1 / 2 * sqrt( (1 + (m == 1)) .* (n + m) .* (n - m + 1) );
        anm(i) = 0;
        bnm = 1 / 2 * sqrt( (n - m) .* (n + m + 1) );
        bnm(i) = sqrt(n0 .* (n0 + 1) / 2);
        SHCs = SHCs .* (n + 2) .* (n + 3);
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = anm .* circshift(pnm, 1) + bnm .* circshift(pnm, -1);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 0 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vyyy'
        i = m == 0; i12 = (m == 1 | m == 2);
        n12 = n(i12); m12 = m(i12);
        anm = 1 / 8 ./ m .* sqrt( (2 * n + 1) ./ (2 * n - 1) .* (1 + (m == 3)) .* (m + n) .* (n + m - 1) .* (n + m - 2) .* (n + m - 3) .* (n - m + 1) .* (n - m + 2) );
        anm(i | i12) = 0;
        bnm = 1 / 8 ./ m .* sqrt( (2 * n + 1) ./ (2 * n - 1) .* (n + m) .* (n + m - 1) ) .* ( 3 * m .* m - m .* (2 * n + 3) + 3 * n .* n + 17 * n + 18 );
        bnm(i) = 0;
        bnm(i12) = 3 ./ (4 * m12 .* (m12 + 1)) .* sqrt( (2 * n12 + 1) ./ (2 * n12 - 1) .* (1 + (m12 == 1)) .* (m12 + n12) .* (n12 + m12 -1) ) .* (n12 + m12 + 1) .* (n12 + m12 + 2);
        cnm = 1 / 8 ./ m .* sqrt( (2 * n + 1) ./ (2 * n - 1) .* (m - n) .* (m - n + 1) ) .* ( 3 * m .* m + m .* (2 * n + 3) + 3 * n .* n + 17 * n + 18 );
        cnm(i) = 0;
        cnm(i12) = 3 ./ (2 * m12 .* (m12 + 2)) .* sqrt( (2 * n12 + 1) ./ (2 * n12 - 1) .* (n12 - m12) .* (n12 - m12 - 1) ) .* (n12 + m12 + 1) .* (n12 + m12 + 2);
        dnm = 1 / 8 ./ m .* sqrt( (2 * n + 1) ./ (2 * n - 1) .* (n - m) .* (n + m + 1) .* (n + m + 2) .* (n - m - 1) .* (n - m - 2) .* (n - m - 3) );
        dnm(i) = 0;
        id = listind - n;
        SHCs = SHCs .* m * -1i;
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = pnm(id);
            pnm = anm .* circshift(pnm, 3) + bnm .* circshift(pnm, 1) + cnm .* circshift(pnm, -1) + dnm .* circshift(pnm, -3);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 1 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vyyz'
        i = m < 2; i0 = m == 0; i1 = m == 1;
        n0 = n(i0); n1 = n(i1);
        anm = -1 / 4 * sqrt( (1 + (m == 2)) .* (n + m) .* (n + m - 1) .* (n - m + 1) .* (n - m + 2) );
        anm(i) = 0;
        bnm = -1 / 2 * ((n + 1) .* (n + 2) + m .* m);
        bnm(i1) = -1 / 4 * (n1 + 2) .* (n1 + 3);
        cnm = -1 / 4 * sqrt( (n - m) .* (n - m - 1) .* (n + m + 1) .* (n + m + 2) );
        cnm(i0) = -1 / 4 * sqrt( 2 * (n0 - 1) .* n0 .* (n0 + 1) .* (n0 + 2) );
        SHCs = SHCs .* -(n + 3);
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = anm .* circshift(pnm, 2) + bnm .* pnm + cnm .* circshift(pnm, -2);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 0 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
    case 'vyzz'
        i = m == 0;
        anm = 1 / 2 ./ m .* sqrt( (2 * n + 1) ./ (2 * n - 1) .* (1 + (m == 1)) .* (n + m) .* (n + m - 1) );
        anm(i) = 0;
        bnm = 1 / 2 ./ m .* sqrt( (2 * n + 1) ./ (2 * n - 1) .* (n - m) .* (n - m - 1) );
        bnm(i) = 0;
        id = listind - n;
        SHCs = SHCs .* -(n + 2) .* (n + 3) .* m * -1i;
        for ii = 1 : nthe
            pnm = Pnm(ii);
            pnm = pnm(id);
            pnm = anm .* circshift(pnm, 1) + bnm .* circshift(pnm, -1);
            pnm = pnm .* SHCs;
            tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
            for mfor = 1 : nmax - 1
                ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
                tp(mfor + 1) = sum(pnm(ind));
            end
            AB(ii, :) = tp;
        end
end

% 2nd step
fld = real(AB) * cosd(mlam) + imag(AB) * sind(mlam);

end
%% subroutine
function u = kapa(n, typ, GM, a, r)

arguments
    n (:, 1)
    typ
    GM (1, 1)
    a (1, 1)
    r (1, 1)
end

ir = 1 / r;
ir2 = ir * ir;
ir3 = ir2 * ir;
ir4 = ir2 * ir2;

logq = log1p((a - r) / r);
qn = exp(n * logq);

switch typ
    case 'shs'
        u = 1;
    case {'v'}
        u = GM * ir;
    case {'vx', 'vy'}
        u = GM * ir2; 
    case 'vz'
        u = -(n + 1) * (GM * ir2);
    case {'vzz'}
        u = (n + 1) .* (n + 2) * (GM * ir3);
    case {'vxx', 'vyy', 'vxy', 'vxz', 'vyz'}
        u = GM * ir3;
    case {'vzzz'}
        u = -(n + 1) .* (n + 2) .* (n + 3) * (GM * ir4);
    case {'vxxx', 'vxxy', 'vxxz', 'vxyy', 'vxyz', 'vxzz', 'vyyy', 'vyyz', 'vyzz'}
        u = GM * ir4;
end

u = u .* qn;

end