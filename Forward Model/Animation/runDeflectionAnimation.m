addpath('..');
clear all, close all, clc

%% Choose plate/wave parameters
alpha = 10;      % non-dim frequency
beta  = 1e-2;    % non-dim rigidity (base chosen as 1e-2)
gamma = 0.1;     % non-dim draught
R     = 1;       % non-dim plate radius
nu    = 0.3;     % Poisson's ratio

% Note: The smaller beta -> more modes of vibration -> increase M,P,N
M     = 10;      % vertical mode truncation
P     = 6;       % Gegenbauer mode truncation
N     = 10;      % angular mode truncation (per ch.2 summary: N~10 at alpha=10)

fprintf('Precomputing amplitudes (Steps 0-6, once)...\n');
data = precomputeDeflectionData(alpha, beta, gamma, R, nu, M, P, N);
fprintf('Done.\n\n');

%% 3D animated surface over one wave period
fprintf('Rendering 3D surface animation...\n');
nr = 60; ntheta = 60; nFrames = 40;
plotDeflectionSurface(data, nr, ntheta, nFrames);

%% 1D time traces at a few chosen points
fprintf('Rendering 1D time-series plot...\n');
points = [0,   0;      % centre
          0.5, 0;      % mid-radius, theta=0
          0.9, 0;      % near edge, theta=0
          0.9, pi];    % near edge, opposite side
nTau = 100;
plotDeflectionTimeSeries(data, points, nTau);