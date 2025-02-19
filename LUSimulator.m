function results = LUSimulator(param)
% LUSimulator simulates the lattice ultrasensitve model of chemotaxis
% signaling. It takes as input a struct (param) containing all the relevant
% parameters for the model. Parameters from fitting are found in
% defaultParams.mat, and custom param structs can be generated using the
% function parameterGeneratorLU.m. It outputs "results", a struct
% containing the all the mean values of methylation, activity, receptor state,
% and control parameter noise at every timestep. It also contains the
% ending lattice configuration of activity and methylation states

% turn off adaptation when param.adapt is set to 0
if param.adapt == 0
    param.R  = 0;
    param.B = 0;
end

if param.kratio == 1
    param.kUp = param.k;
    param.kDown = param.k;
end

% set starting parameters
A = param.Astart;
M = param.Mstart;
cExp = param.cExp0;

N = length(param.Astart); %lattice size
sweepfrac = 0.1; %maximum fraction of active/inactive kinases that can change state in one sweep

convMat = [0,1,0;1,0,1;0,1,0]; % used with convolution to calculate sum of neighboring activities

% determine maximum reasonable number of timesteps in simulation, for preallocation
if param.adapt == 1
    zMax = (1+exp(-param.alphaS*param.m0r+4*param.alphaS)*(1 + param.kUp) + param.R + param.B)*N^2; 
else
    zMax = (1+exp(-param.alphaS*param.m0r+param.alphaS*mean(param.Mstart,'all'))*(1 + param.kUp))*N^2; 
end
dtMin = 1/zMax;
lenMax = ceil(param.tf/dtMin/500); % estimate of maximum length, reduced to a more reasonable size

% preallocate variables
results.ts   = zeros(lenMax,1);
results.As   = zeros(lenMax,1);
results.Ms   = zeros(lenMax,1);
results.ps   = zeros(lenMax,1);
results.cExp = zeros(lenMax,1);
% results.actCluster = zeros(lenMax,5);
% results.inactCluster = zeros(lenMax,5);
results.actCluster = zeros(lenMax,1);

zMat = cell(4,1);
zLin = zMat;

t = 0; % start at time 0
counter = 0; % track number of time steps
while t<param.tf % continue until simulation has run for a time param.tf
    counter = counter + 1; % increment number of timesteps
    
    % calculate mean receptor state
    % f = param.Ea0-param.alpha*M + log((1+param.L/param.Ki)/(1+param.L/param.Ka));
    f = param.receptorHillCoef*(-param.alphaR*(M - param.m0r) + log((1+param.L/param.Ki)/(1+param.L/param.Ka)));
    p = 1./(1+exp(f));
    
    g = param.alphaS*(M-param.m0s); %calculate spreading modification

    % determine number of active and inactive kinases adjacent to each
    % lattice site
    if param.continuousBoundary == 1
        % expand A for continuity
        continuousA = [A A A;A A A;A A A];
        continuousA = continuousA(N:2*N+1,N:2*N+1);
        adjacenciesAct   = conv2(continuousA,convMat,'same')/4;
        adjacenciesInact = conv2(1-continuousA,convMat,'same')/4;
        adjacenciesAct   = adjacenciesAct(2:end-1,2:end-1);
        adjacenciesInact = adjacenciesInact(2:end-1,2:end-1);
    else
        adjacenciesAct   = conv2(A,convMat,'same')/4;
        adjacenciesInact = conv2(1-A,convMat,'same')/4;
    end
    
    % record results
    Atot = sum(A,'all');
    Aavg = Atot/N^2;
    Mavg = sum(M,'all')/N^2;
    pavg = sum(p,'all')/N^2;
    
    results.As(counter) = Aavg;
    results.ps(counter) = pavg;
    results.Ms(counter) = Mavg;
    results.ts(counter) = t;  
    results.cExp(counter) = cExp;
    % for i = 1:5
    %     results.actCluster(counter,i) = sum(adjacenciesAct == (i-1)/4,'all')/N^2;
    %     results.inactCluster(counter,i) = sum(adjacenciesAct == (i-1)/4,'all')/N^2;
    % end
    results.actCluster(counter) = sum(adjacenciesAct.*(1-A),'all')/N^2;

    % Calculate Propensity Functions for the following reactions at each site
    % 1: increase A
    % 2: decrease A
    % 3: increase M
    % 4: decrease M
    
    if param.definedc == 1 % use if running sim at a set value of c
        % The choice function of c in front of these terms is equivalent to
        % assuming p = 0.5 for all values of c ~= 1.
        zMat{1} = (param.c+1)/2*(adjacenciesAct + param.kUp).*(1-A);
        zMat{2} = (1/param.c+1)/2*(adjacenciesInact + param.kDown).*(A);
    elseif param.frankModification == 1 % used when spontaneous switching is independent of methylation
        zMat{1} = p.*(adjacenciesAct.*(1+exp(g)) + param.kUp).*(1-A);
        zMat{2} = (1-p).*(adjacenciesInact.*(1+exp(-g)) + param.kDown).*(A);
    else % used when input is based on ligand and methylation
        zMat{1} = p.*(1+exp(g)).*(adjacenciesAct + param.kUp).*(1-A);
        zMat{2} = (1-p).*(1+exp(-g)).*(adjacenciesInact + param.kDown).*(A);
    end
    zMat{1} = exp(cExp)*zMat{1}; % include the impact of external noise in the control parameter
    zMat{3} = param.R*(1-A);
    zMat{4} = param.B*A*(1 + (Aavg-1)*param.globalAdapt);

    % transform matrices of propensities into a linear sum which can be
    % used to weight the random selection of the next reaction
    for i = 1:length(zMat)
        zLin{i} = reshape(zMat{i},[N^2,1]);
    end
    zLinAll = [zLin{1};zLin{2};zLin{3};zLin{4}];
    zLinSum = cumsum(zLinAll);
    z = zLinSum(end);

    % determine number of reactions per time step
    if param.tauLeap == 1 % if using tau leaping
        numPerSweep = max(floor(sweepfrac*min(Atot,N^2-Atot)),1);
    else % if not using tau leaping
        numPerSweep = 1;
    end

    randNumsReact  = rand(1,numPerSweep); % pregen random numbers
    
    dt = numPerSweep/z*log(1/rand()); % determine timestep size
    t = t+dt; % increment time

    cExp = cExp - param.externalTimeScale*cExp*dt + param.externalNoiseScale*(randn)*sqrt(dt); % update control parameter noise

    for i = 1:numPerSweep
        index = find(zLinSum>randNumsReact(i)*z,1); % determine type and location of next reaction
        if index <= N^2 % make kinase active
            [x,y] = ind2sub([N,N],index);
            A(x,y) = 1;
        elseif index <= 2*N^2 % make kinase inactive
            [x,y] = ind2sub([N,N],index-N^2);
            A(x,y) = 0;
        elseif index <= 3*N^2 % increase methylation if methylation is below the maximum
            [x,y] = ind2sub([N,N],index-2*N^2);
            if M(x,y) < 4
                M(x,y) = M(x,y)+param.dM;
            end
        elseif index <= 4*N^2 % decrease methylation if methylation is above the minimum
            [x,y] = ind2sub([N,N],index-3*N^2);
            if M(x,y) > 0
                M(x,y) = M(x,y)-param.dM;
            end
        end
    end

end

% record final lattice state
results.Aend = A;
results.Mend = M;

% remove unused preallocations
results.ts = results.ts(1:find(results.ts,1,'last'));
results.As = results.As(1:length(results.ts));
results.Ms = results.Ms(1:length(results.ts));
results.ps = results.ps(1:length(results.ts));
results.cExp = results.cExp(1:length(results.ts));
results.actCluster = results.actCluster(1:length(results.ts));

end