function C00 = shell_C00_cons(r1, r2, rho, M)
%{
    caculate a potential of a spherical shell with constant density C00
Input:
    r1: inner radius [m]
    r2: outer radius [m]
    rho: density [kg/m^3]
    M: scaling constant of mass [m^3 kg^-1 s^-2]
Return:
    C00: (1, 1) the C00 of SHCs
%}
arguments
    r1 (1, 1)
    r2 (1, 1)
    rho (1, 1)
    M (1, 1)
end

digits(100);
C00 = vpa(4 / 3 * vpa(pi, 32)  * rho * (r2 ^ 3 - r1 ^ 3) / M);
C00 = double(C00);

end
