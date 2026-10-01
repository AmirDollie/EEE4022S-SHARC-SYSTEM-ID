function P = setupStudyPaths()
%SETUPSTUDYPATHS  Add every folder a Design-and-Robustness study needs to the path; return key folders.
%   P = setupStudyPaths()
%   P.root (repo), P.dr (Design and Robustness), P.results (Parameter Inversion/Results),
%   P.cache (Results/twinCache)
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.
P.dr = fileparts(mfilename('fullpath'));
piRoot = fileparts(P.dr);
P.root = fileparts(fileparts(piRoot));
addpath(P.dr);
addpath(fullfile(P.dr, 'Engine'));
addpath(fullfile(P.dr, 'Reporting'));
addpath(fullfile(P.root, 'Transmissibility'));
addpath(fullfile(P.root, 'Inverse Mapping'));
addpath(fullfile(P.root, 'Forward Model'));
addpath(fullfile(P.root, 'Forward Model', 'Animation'));
addpath(fullfile(P.root, 'SSI'));
P.results = fullfile(piRoot, 'Results');
P.cache = fullfile(P.results, 'twinCache');
P.imuResults = fullfile(P.root, 'IMU Characterisation', 'Results');
end