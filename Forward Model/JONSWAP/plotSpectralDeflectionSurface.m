function plotSpectralDeflectionSurface(specData, nr, ntheta, nFrames, tMax, exaggeration)
%PLOTSPECTRALDEFLECTIONSURFACE Animated 3D surface of JONSWAP-forced plate
%deflection. x,y plotted in non-dimensional r,theta (matching the
%original single-frequency plotDeflectionSurface.m style). Colour always
%reflects the TRUE physical deflection in mm; the plotted HEIGHT is
%artificially exaggerated (standard practice for genuinely small
%deflections relative to plate radius) and the exaggeration factor is
%printed clearly in the title so it is never mistaken for true scale.
%
%   PLOTSPECTRALDEFLECTIONSURFACE(SPECDATA, NR, NTHETA, NFRAMES, TMAX, EXAGGERATION)
%   EXAGGERATION - optional; if omitted, auto-chosen so the tallest
%       feature occupies about 30% of the plate radius on screen.

    thisDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(thisDir, '..'));
    addpath(fullfile(thisDir, '..', 'Animation'));

    R = specData.R;   % non-dim, for x,y

    r = linspace(0, R, nr);
    theta = linspace(0, 2*pi, ntheta);
    [rGrid, thetaGrid] = meshgrid(r, theta);
    xGrid = rGrid .* cos(thetaGrid);   % non-dim (matches original style)
    yGrid = rGrid .* sin(thetaGrid);

    numBins = length(specData.binData);
    etaGridAll = zeros(size(rGrid,1), size(rGrid,2), numBins);
    fprintf('Precomputing spatial pattern for %d bins...\n', numBins);
    for i = 1:numBins
        etaGridAll(:,:,i) = evaluateDeflection(specData.binData{i}, rGrid, thetaGrid);
    end

    tVals = linspace(0, tMax, nFrames);
    zetaAll_m = zeros([size(rGrid), nFrames]);   % TRUE deflection, metres
    for i = 1:numBins
        phase = exp(1i*(specData.omega(i)*tVals + specData.epsilon(i)));
        zetaAll_m = zetaAll_m + specData.a(i) * real(etaGridAll(:,:,i) .* reshape(phase, 1, 1, nFrames));
    end
    zetaAll_mm = zetaAll_m * 1000;   % for the colorbar, in mm

    if nargin < 6 || isempty(exaggeration)
        maxAbsZeta = max(abs(zetaAll_m(:)));
        exaggeration = (0.3 * R) / maxAbsZeta;
    end
    plottedHeight = zetaAll_m * exaggeration;
    heightLimit = max(abs(plottedHeight(:)));
    colorLimit_mm = max(abs(zetaAll_mm(:)));

    figure;
    for f = 1:nFrames
        surf(xGrid, yGrid, plottedHeight(:,:,f), zetaAll_mm(:,:,f), 'EdgeColor', 'none');
        colormap(divergingColormap);
        clim([-colorLimit_mm, colorLimit_mm]);
        zlim([-heightLimit, heightLimit]);
        cb = colorbar; cb.Label.String = '\zeta (mm, true scale)';
        xlabel('r/H (non-dim)'); ylabel('r/H (non-dim)');
        zlabel('plotted height (exaggerated)');
        title(sprintf('Spectral deflection, t=%.2fs  [vertical exaggeration: %.0fx]', ...
            tVals(f), exaggeration));
        axis equal;
        drawnow;
        pause(0.05);
    end
end

% ------------------------------------------------------------------
function cmap = divergingColormap(n)
    if nargin < 1, n = 256; end
    half = floor(n/2);
    blueToWhite = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
    whiteToRed  = [ones(n-half,1), linspace(1,0,n-half)', linspace(1,0,n-half)'];
    cmap = [blueToWhite; whiteToRed];
end