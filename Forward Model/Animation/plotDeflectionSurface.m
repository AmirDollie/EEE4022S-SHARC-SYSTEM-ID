function plotDeflectionSurface(data, nr, ntheta, nFrames)
%PLOTDEFLECTIONSURFACE Animated 3D surface of plate deflection over one
%full non-dimensional wave period (tau in [0, 2*pi)).
    addpath('..');

    R = data.R;
    r = linspace(0, R, nr);
    theta = linspace(0, 2*pi, ntheta);
    [rGrid, thetaGrid] = meshgrid(r, theta);
    xGrid = rGrid .* cos(thetaGrid);
    yGrid = rGrid .* sin(thetaGrid);

    etaStatic = evaluateDeflection(data, rGrid, thetaGrid);
    zLimit = max(abs(etaStatic(:)));

    tauVals = linspace(0, 2*pi, nFrames+1);
    tauVals(end) = [];

    figure;
    for f = 1:length(tauVals)
        zGrid = real(etaStatic .* exp(-1i*tauVals(f)));
        surf(xGrid, yGrid, zGrid, 'EdgeColor', 'none');
        colormap(parula); colorbar;
        zlim([-zLimit, zLimit]);
        xlabel('x'); ylabel('y'); zlabel('\eta (non-dim)');
        title(sprintf('Plate deflection, \\tau = %.2f', tauVals(f)));
        axis equal;
        drawnow;
        pause(0.05);
    end
end