function C00 = shell_C00_poly(r1, r2, rho, M, r0)
%{
    caculate a potential of a spherical shell with radial polynomial density C00
Input:
    r1: inner radius [m]
    r2: outer radius [m]
    rho: density [kg/m^3] [p0, p1, p2,...,pP], Polynomial coefficient
    M: M: scaling constant of mass [m^3 kg^-1 s^-2]
Return:
    C00: (1, 1) the C00 of SHCs
%}
arguments
    r1 (1, 1)
    r2 (1, 1)
    rho (1, :)
    M (1, 1)
    r0 (1, 1)
end
digits(100);
r1 = vpa(r1);
r2 = vpa(r2);
r0 = vpa(r0);

P = length(rho) - 1;
C00 = 0;
t2 = r2 - r0; t1 = r1 - r0; 
for p = 0 : P
    t = (t2 ^ (p + 3) - t1 ^ (p + 3)) / (p + 3) + 2 * r0 * (t2 ^ (p + 2) - t1 ^ (p + 2)) / (p + 2) + r0 ^ 2 * (t2 ^ (p + 1) - t1 ^ (p + 1)) / (p + 1);
    C00 = C00 + rho(p + 1) * t;
end
C00 = 4 * pi / M * C00;
C00 = double(C00);
end
