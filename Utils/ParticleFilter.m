% Standard multiple or single model particle filter. Assumes nearly
% constant velocity model for state transition and range measurements.
%
% Input:
%   particles  - A structure consisting of states, weights, and models of
%                each particle.
%   sensors    - Locations of the anchor nodes
%   meas       - Range measurement to each anchor.
%   Parameters - A structure consisting of various parameters.
%
% Output:
%   particles  - Updated particles.
%   state_estimate - Mean of the updated state.
%
%

function [particles, state_estimate] = ParticleFilter(particles, sensors, meas, Parameters)

nParticles = Parameters.PF_NumberOfParticles;
statedim = 2*Parameters.SpaceDimension;

if strcmp(Parameters.PF_LikelihoodFunction,'Gaussian')
    Rk = Parameters.VarianceOfRangeMeasurement * eye(length(meas));
    normfac = 1/sqrt(det(2*pi*Rk)); % Gaussian normalizing factor
elseif strcmp(Parameters.PF_LikelihoodFunction,'Top-Hat')
    % likelihood with top hat distribution
    MIN_ERR = Parameters.PF_Likel_MinimumError;
    MAX_ERR = Parameters.PF_Likel_MaximumError;
end

if (length(unique(particles.models)) > 1)
    % multiple model version, propose new models first
    pij_c = cumsum(Parameters.PF_ModelTransitionProbMat, 2);
    for k1 = 1:nParticles
        
        cMode = particles.models(k1);
        u = rand;
        modProb = pij_c(cMode,:);
        idx = u>modProb;
        idx = find(idx==0,1);
        
        particles.models(k1) = idx;
    end
end

% now propose new state vectors - prediction, conditioned on the mode
for k1 = 1:nParticles
    q1 = Parameters.ProcessNoiseIntensity(particles.models(k1));
    [Phi,Qk] = GetMotionModel(Parameters.SpaceDimension,'StateTransitionMatrix_CV',Parameters.SamplingTime,'ProcessNoiseCovarianceMatrix_CV',Parameters.SamplingTime,q1);
    particles.states(1:statedim,k1) = Phi * particles.states(1:statedim,k1) + chol(Qk) * randn(statedim,1);
end

if ~isempty(meas)    
    % update the weights
    if (Parameters.SpaceDimension == 2)
        sPtsPos = [particles.states(1,:);particles.states(3,:)];
    elseif (Parameters.SpaceDimension == 3)
        sPtsPos = [particles.states(1,:);particles.states(3,:);particles.states(5,:)];
    end
    
    numSens = length(meas);
    
    for k1 = 1:nParticles
        
        xydiff = repmat(sPtsPos(:,k1)',numSens,1) - sensors;
        predM  = sqrt(sum(xydiff.^2,2));    
        nuk = predM - meas;
        
        if strcmp(Parameters.PF_LikelihoodFunction,'Gaussian')
            % likelihood with Gaussian distribution
            particles.weights(:,k1) = normfac * exp(-0.5 * nuk' * (Rk \ nuk) ) + 1e-99;
        elseif strcmp(Parameters.PF_LikelihoodFunction,'Top-Hat')
            error_grid = nuk;
%             cweight = 0.01^sum(((error_grid > MAX_ERR) | (error_grid < MIN_ERR)));
%             cweight = cweight * prod((error_grid(((error_grid <= 0) & (error_grid > MIN_ERR))) - MIN_ERR)*(1/MIN_ERR) + 0.01);
%             cweight = cweight * prod((- error_grid(((error_grid > 0) & (error_grid <= MAX_ERR))) + MAX_ERR)*(1/MAX_ERR) + 0.01);
%             particles.weights(1,k1) = cweight;
            particles.weights(1,k1) = prod(((error_grid < MAX_ERR) & (error_grid > MIN_ERR)) + 0.01);
        else
            error('Unknown likelihood distribution.');
        end
    end
    % normalize weights
    particles.weights = particles.weights/sum(particles.weights);
    
%     % find current state estimate
%     max_weight = max(particles.weights);
%     prob2 = (particles.weights > (1e-4 * max_weight)) .* particles.weights;
%     summed_probability = sum(prob2);
%     prob2 = prob2(ones(4,1),:);
%     state_estimate = sum(particles.states .* prob2,2)/summed_probability;
    
    % resample
    particlesToKeep = ResamplingAlgorithms(particles.weights,'SystematicResampling');    
    
    particles.weights = (1/nParticles)*ones(1,nParticles);
    particles.states = particles.states(:,particlesToKeep);
    particles.models = particles.models(:,particlesToKeep);
    
    state_estimate = mean(particles.states,2);    
else
    state_estimate = mean(particles.states,2);
end


