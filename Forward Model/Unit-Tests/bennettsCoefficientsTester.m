% Tester to see if the Bennetts reduction map actually converges to the
% direct method for phi
clear all, close all, clc

%alpha, beta, gamma test cases
cases = {
    10, 1e-2, 0.1;
    10, 1e-4, 1e-3;
    5,  1e-3, 0.05;
};


zTest = [-0.5, -0.9];   % interior points


for c = 1:size(cases,1)
    alpha = cases{c,1}; beta = cases{c,2}; gamma = cases{c,3};
    fprintf('\n=== Case %d: alpha=%.4g, beta=%.4g, gamma=%.4g ===\n', c, alpha, beta, gamma);
    
    %arb set of lengths. We need to show the difference between the
    %reconstruction and the original decrease as M increases
    for M = [8, 20, 40, 80]
        %calc the disp roots
        kappa = dispersionRoots(alpha, beta, gamma, M);
        %construct V
        V = bennettsCoefficients(alpha, beta, gamma, kappa);
        %separate the m>0 cases
        kappaTower = kappa(3:end);

        for j = 1:2
            kappaJ = kappa(j);   % V(:,1) <-> kappa(1)=kappa_{-2}; V(:,2) <-> kappa(2)=kappa_{-1};

            % build the full list of z-points to check: the edge, plus
            % every interior point in zTest
            zPoints = [-gamma - 1e-6, zTest];

            for zi = 1:length(zPoints)
                zVal = zPoints(zi);

                direct = verticalEigenfunction(zVal, kappaJ, gamma);
                reconstructed = 0;
                for m = 1:length(kappaTower)
                    reconstructed = reconstructed + V(m,j)*verticalEigenfunction(zVal, kappaTower(m), gamma);
                end

                fprintf('  M=%3d, j=%d, z=%6.3f: direct=%.6f%+.6fi, reconstructed=%.6f%+.6fi, diff=%.3e\n', ...
                    M, j, zVal, real(direct), imag(direct), real(reconstructed), imag(reconstructed), ...
                    abs(direct-reconstructed));
            end
        end
    end
end