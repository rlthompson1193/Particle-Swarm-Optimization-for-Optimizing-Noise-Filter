clear all;

%~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters for adjustment
%~~~~~~~~~~~~~~~~~~~~~~~~~~
np = 1000; % Number of particles
order = 4; % Filter order
UB = 10; % Filter coefficient upper bound
LB = -10; % Filter coefficient lower bound
iterations = 750; % Number of iterations
inertialCoefficient = 0.008; % Inertial scale factor
pBestInfluence = 0.008; % Personal best scale factor
gBestInfluence = 0.01; % Global best scale factor
% Time vector, default set for 100 points on the interval [0 1]
t = linspace(0,1,100);

% Create discrete ideal input signal vector
% For simplicity, this is a simple 2 Hz sine wave
% with an amplitude of 1
cleanSignal = sin(2*pi*2*t);
%plot(t,cleanSignal)
% Create noisy version of input signal vector
% This can be customized. Currently loading a vector
% of noise which is uniformly distributed on [-1 1]
% then scaling so noise amplitude interval is 10%
% of signal amplitude
noiseVector = csvread('noiseVector.csv');
noisySignal = cleanSignal + 0.1 * noiseVector;
%plot(t,noisySignal)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% After this line, no code should need to be altered. All
% remaining code is parameterized based on the above variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initialize random starting points for particles
% Each column will contain the coefficients for each particle
CurrentParticlePositions = rand(order + 1, np) * (UB - LB) +…
 (LB * ones(order + 1, np));
NewParticlePositions = CurrentParticlePositions;
PreviousParticlePositions = CurrentParticlePositions;

% Run filter on initial values, then initialize personal best and global
% best vectors
currentASE = zeros(1,np);
newASE = currentASE;
for i = 1:np
 currentASE(i) = ComputeASE(cleanSignal, noisySignal,
NewParticlePositions(:,i));
end
pBest = currentASE;
pBestPositions = NewParticlePositions;
[gBest, minIndex] = min(pBest);
gBestPositions = NewParticlePositions(:,minIndex);




% Main algorithm loop
for i = 1:(iterations - 1)

 % Find contribution from personal best
 pBestContribution = pBestPositions - CurrentParticlePositions;

 % Find contribution from global best
 for j = 1:np
 gBestContribution(:,j) = gBestPositions -
CurrentParticlePositions(:,j);
 end

 % Find contribution from intertia
 inertialContribution = CurrentParticlePositions -
PreviousParticlePositions;

 % Find next particle positions, but don't move particles yet
 NewParticlePositions = CurrentParticlePositions + (pBestContribution *
pBestInfluence) ...
 + (gBestContribution * gBestInfluence) + ...
 (inertialContribution * inertialCoefficient);

 % Remember current positions for future inertia calculation
 PreviousParticlePositions = CurrentParticlePositions;

 % Move particles to new positions
 CurrentParticlePositions = NewParticlePositions;

 % Find ASE for each current particle position and check it against global
 % and personal bests, then update those if appropriate
 for j = 1:np
 % Calculate ASE for particle j
 ASEj = ComputeASE(cleanSignal, noisySignal,
CurrentParticlePositions(:,j));
 % Update personal best if appropriate
 if (ASEj < pBest(j))
 pBest(j) = ASEj;
 pBestPositions(:,j) = CurrentParticlePositions(:,j);
 end
 % Update global best if appropriate
 if (ASEj < gBest)
 gBest = ASEj;
 gBestPositions = CurrentParticlePositions(:,j);
 end
 end
end
% display best solution found by algorithm
gBest
transpose(gBestPositions)
%pBestPositions
% This function passes the noisy signal through the desired filter,
% compares the output to the original clean signal, and computes the ASE
function ASE = ComputeASE(cleanSignal, noisySignal, coefficients)

 iterations = length(noisySignal) - length(coefficients) + 1;
 outputSignal = zeros(1, iterations);

 % Calculate and build output signal
 for i=1:iterations
 outputSignal(i) = flip(noisySignal(i:(i + length(coefficients) - 1)))
...
 * coefficients;
 end

 % Compute error
 errorVector = cleanSignal(length(coefficients):end) - outputSignal;
 ASE = sum(errorVector .^ 2) / length(errorVector);

end
