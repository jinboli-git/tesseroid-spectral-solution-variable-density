function [V, Vz, Vzz, Vxx, Vzzz, Vxxz] = shell_geff_cons(R1, R2, r, rho, G)
%{
    caculate a gravity field of a spherical shell with constant density
Input:
    R1: inner radius [m]
    R2: outer radius [m]
    r: the distance from the center at the point [m]
    rho: density [kg/m^3]
    G: gravitational constant
Return:
    V, Vz, Vzz, Vxx, Vzzz, Vxxz
Reference:
    Deng X-L (2023) Evaluation of gravitational curvatures for a tesseroid and spherical shell with arbitrary-order polynomial density. J Geod 97:18. https://doi.org/10.1007/s00190-023-01708-2
%}
arguments
    R1 (1, 1)
    R2 (1, 1)
    r (:, 1)
    rho (1, 1)
    G (1, 1)
end

digits(100);
M0 = vpa(4 / 3 * pi * G * rho * (R2 ^ 3 - R1 ^ 3));

V = M0 ./ r;
Vz = -V ./ r;
Vzz = -2 * Vz ./ r;
Vxx = Vzz / -2; 
Vzzz = -3 * Vzz / r;
Vxxz = Vzzz / -2;

V = double(V);
Vz = double(Vz);
Vzz = double(Vzz);
Vxx = double(Vxx);
Vzzz = double(Vzzz);
Vxxz = double(Vxxz);

end

