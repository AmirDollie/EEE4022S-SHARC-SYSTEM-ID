function plotDeflectionTimeSeries(data, points, nTau)
%PLOTDEFLECTIONTIMESERIES Plots eta(r,theta,tau) vs tau at a few chosen
%(r,theta) points over one full non-dim period.
    addpath('..');

    tauVals = linspace(0, 2*pi, nTau);
    figure; hold on;
    legendEntries = cell(size(points,1),1);

    for p = 1:size(points,1)
        r = points(p,1); theta = points(p,2);
        etaStatic = evaluateDeflection(data, r, theta);
        etaTime = real(etaStatic * exp(-1i*tauVals));
        plot(tauVals, etaTime, 'LineWidth', 1.5);
        legendEntries{p} = sprintf('r=%.2f, \\theta=%.2f', r, theta);
    end

    xlabel('\tau (non-dim phase)'); ylabel('\eta (non-dim)');
    title('Plate deflection vs phase at selected points');
    legend(legendEntries);
    grid on;
end