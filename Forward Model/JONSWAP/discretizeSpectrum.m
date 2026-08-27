function [omega, a, epsilon] = discretizeSpectrum(Hs, w0, gammaJONSWAP, wMin, wMax, numBins)

    %Transforms the continuous JONSWAP spectrum into a finite set of
    %frequency components (amplitude, random-phase)
    

    % Hs = significant wave height (m)
    %change the SWH when necessary!
    
    % w0 = peak frequency (rad/s)
    % gammaJONSWAP = JONSWAP peakedness factor (usually = 3.3)
    % wMin, wMax = frequency range to discretize over. Need wMin > 0
    % numBins = number of discrete frequency components

    %OUTPUTS which are all length(numBins) x 1 vectors:
    % omega = bin-centre angular frequs (rad/s)
    % A = each component's amplitude (from EQUAL VARIANCE MATCHING)
    % epsilon = independent random phases, uniform over [0, 2pi)
    

    % Requires the MSS toolbox's waveSpectrum.m function

    domega = (wMax - wMin) / numBins;
    omega = wMin + ((1:numBins)' - 0.5) * domega;   % bin midpoints

    S = waveSpectrum(7, [Hs, w0, gammaJONSWAP], omega, 0);

    a = sqrt(2 * S * domega);
    epsilon = rand(numBins, 1) * 2 * pi;

end