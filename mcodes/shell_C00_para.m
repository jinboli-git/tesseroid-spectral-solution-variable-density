function C00 = shell_C00_para(r1, r2, rho, M, r0)
%{
    caculate a potential of a spherical shell with radial parabolic density C00
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
r0 = vpa(r0);

mu = rho(2); rho = rho(1);
t2 = rho - mu * (r2 - r0); t1 = rho - mu * (r1 - r0); t0 = rho + mu * r0;

C00 = -t0 ^ 2 * (t2 ^ (-1) - t1 ^ (-1)) + t2 - t1 -2 * t0 * log(t2 / t1);
C00 = 4 * pi / M  * C00 * -(rho / mu) ^ 3;

C00 = double(C00);

end
