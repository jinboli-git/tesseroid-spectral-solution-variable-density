function C00 = shell_C00_exp(r1, r2, rho, M, r0)
%{
    caculate a potential of a spherical shell with radial exponential density C00
Input:
    r1: inner radius [m]
    r2: outer radius [m]
    rho: density [kg/m^3] [rho, mu]
    M: M: scaling constant of mass [m^3 kg^-1 s^-2]
Return:
    C00: (1, 1) the C00 of SHCs
%}
arguments
    r1 (1, 1)
    r2 (1, 1)
    rho (1, 2)
    M (1, 1)
    r0 = 0
end
digits(100);
r1 = vpa(r1);
r2 = vpa(r2);
rho = vpa(rho);
mu = rho(2);
e1 = exp(-mu * (r1 - r0)); e2 = exp(-mu * (r2 - r0));
C00 = (e1 * r1 ^ 2 - e2 * r2 ^ 2) / mu + 2 / (mu ^ 2) * (e1 * r1 - e2 * r2) ...
+ 2 / (mu ^ 3) * (e1 - e2);
C00 = 4 * pi / M * rho(1) * C00;

C00 = double(C00);
end
