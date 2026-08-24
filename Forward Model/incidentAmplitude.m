function AIn = incidentAmplitude(alpha, n, k0, R, M)
%INCIDENTAMPLITUDE The incident-wave amplitude vector A^(In)_n (eq. 2.22).
%   Only the m=0 (travelling) tower mode is forced; every other tower
%   mode gets zero, since the incident wave's own vertical profile
%   matches phi^(0)_0(z) exactly.
    AIn = zeros(M+1, 1);
    AIn(1) = besseli(n, k0*R) / (1i*alpha);
end