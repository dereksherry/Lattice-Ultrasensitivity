function param = parameterGeneratorLU(varargin)
if mod(length(varargin),2)~=0
    error('Incorrect Input Values')
end

for i = 1:2:length(varargin)
    if isequal(varargin{i},'size')
        n = varargin{i+1};
    elseif isequal(varargin{i}, 'nu')
        param.nu = varargin{i+1};
    elseif isequal(varargin{i}, 'kratio')
        param.kratio = varargin{i+1};
    end
end
if ~exist('n','var')
    n = 2^4;
end
if ~exist('param','var') || ~isfield(param, 'nu')
    param.nu = 100; %t(s) = t(simTime)/nu
end
if ~exist('param','var') || ~isfield(param, 'kratio')
    param.kratio = 1;
end

param.k  = 0.2;
param.kUp = 2*param.k*param.kratio/(1+param.kratio);
param.kDown = 2*param.k*1/(1+param.kratio);
param.c = ((1+2*param.kDown)/(1+2*param.kUp));

param.adapt = 0; % include adaptation?
param.tauLeap = 1; % perform tau Leaping?
param.globalAdapt = 1; % include an additional factor of average activity in demethylation
param.homogeneous = 0; % determines whether all methylation changes in sync (this function has not been updated to the current model
param.definedc = 0; % requires externalfluctuations to be on off function
param.continuousBoundary = 0; % 1 for periodic boundaries
param.forScaling = 0;
param.frankModification = 1;

param.externalTimeScale  = 10^-3/param.nu;
externalNoiseVariance = 0.1;
externalNoiseVariance = 0;
param.externalNoiseScale = sqrt(externalNoiseVariance*2*param.externalTimeScale);

param.receptorHillCoef = 1;
param.receptorHillCoef = 1.75;
param.m0r    = -2; % receptor 0 methylation effect
param.m0r    = 1.7/2;
param.alphaR = 1; % Receptor methylation dependence
param.alphaR = 0.9*2;
param.m0s    = -1; % spreading 0 methylation effect
param.m0s    = -6.7/2;
param.alphaS   = 1; % spreading methylation dependence
param.alphaS   = 0.1*2;
param.dM    = 1; % how much the methylation value of a receptor changes by in one reaction.
param.dM    = 1/12;


param.totalAdapt  = 0.02/param.nu/param.dM; %methylation rate
param.RBratio     = 0.8;
A0 = 0.3;
param.RBratio = A0^2/(1-A0);
param.R = param.RBratio/(param.RBratio+1) * param.totalAdapt;
param.B = 1/(param.RBratio+1) * param.totalAdapt;

param.Ki    = 0.0062; %uM
param.Ka    = 1;

%From Shimizu Fit
param.Ki    = 4*10^-5;
param.Ka = 0.73;

param.Ka = 0.5;
param.Ki = 0.01;

% param.L = 0.074;
% param.L = 100;
param.L = 0.0;

param.cExp0 = 0;
param.Astart = randi(2,n)-1;
param.Mstart = randi(2,n)+1;
param.Mstart = 1*ones(n);
if param.homogeneous == 1
    param.Mstart = 2*ones(n);
end
%param.M0 = 2*ones(n);

param.tf = 10*param.nu;