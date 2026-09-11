function y = syntheticSensorData(specData, sensorLocations, tVec)
%SYNTHETICSENSORDATA Builds a multi-channel time series matrix from
%precomputed Forward Model deflection data, for feeding into
%buildHankelMatrix.m and the rest of the SSI pipeline.
%
%   Y = SYNTHETICSENSORDATA(SPECDATA, SENSORLOCATIONS, TVEC)
%
%   SPECDATA        - struct from precomputeSpectralData.m
%   SENSORLOCATIONS - Nsensors x 2 matrix, each row [r, theta] for one
%                     sensor location (non-dimensional r, must be <= R)
%   TVEC            - 1 x n vector of real physical times (seconds)
%
%   Returns Y, an Nsensors x n matrix: row i is
%   evaluateSpectralDeflection(specData, sensorLocations(i,1),
%   sensorLocations(i,2), tVec), i.e. exactly the shape
%   buildHankelMatrix.m expects (l channels x n samples).
%                 
%   This is pure orchestration: it calls the already-validated
%   evaluateSpectralDeflection.m once per sensor and stacks the results
%   as rows!

    thisDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(thisDir, '..', 'Forward Model'));
    addpath(fullfile(thisDir, '..', 'Forward Model', 'Animation'));
    addpath(fullfile(thisDir, '..', 'Forward Model', 'JONSWAP'));

    nSensors = size(sensorLocations, 1);
    y = zeros(nSensors, length(tVec));

    for i = 1:nSensors
        r = sensorLocations(i,1);
        theta = sensorLocations(i,2);
        y(i,:) = evaluateSpectralDeflection(specData, r, theta, tVec);
    end
end