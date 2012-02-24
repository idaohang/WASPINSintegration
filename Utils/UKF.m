% An implementation of the unscented Kalman filter 
% Input:
%   Xkk        - Prior state.
%   Pkk        - Prior state covariance.
%   meas       - Range measurements to anchors.
%   sensors    - Locations of the anchors.
%   Parameters - A structure consisting of various parameters.
%
% Output:
%   Xk1        - Updated state.
%   Pk1        - Updated state covariance
%   likeli     - Filter calculated likelihood.
%
%

function [Xk1,Pk1,likeli] = UKF(Xkk, Pkk, meas, sensors, Parameters)

L = length(Xkk);
alpha = 5e-1;
kappa = 0;
lambda = alpha^2*(L + kappa) - L;
mybeta = 2;

num_sigmaPts = 2*L + 1;

[Phi,Q] = GetMotionModel(Parameters.SpaceDimension,'StateTransitionMatrix_CV',Parameters.SamplingTime, ...
    'ProcessNoiseCovarianceMatrix_CV',Parameters.SamplingTime,Parameters.ProcessNoiseIntensity);

xkminus = Phi * Xkk;
xkminusmat = repmat(xkminus,1,2*L+1);
Pkminus = Phi * Pkk * Phi' + Q;

if isempty(meas)
    % no measurements to update
    Xk1 = xkminus;
    Pk1 = Pkminus;
    likeli = 1;
else
    % have measurements proceed
    sqrtPkk = chol((L+lambda)*Pkminus)';
    
    sigmaPts = [zeros(size(xkminus)) sqrtPkk -sqrtPkk];
    sigmaPts = xkminusmat + sigmaPts;
    
    if (Parameters.SpaceDimension == 2)
        sPtsPos = [sigmaPts(1,:);sigmaPts(3,:)];
    elseif (Parameters.SpaceDimension == 3)
        sPtsPos = [sigmaPts(1,:);sigmaPts(3,:);sigmaPts(5,:)];
    end
    numSens = size(sensors,1);
    Yk1k = zeros(numSens,num_sigmaPts);
    
    for k1=1:num_sigmaPts
        xydiff = repmat(sPtsPos(:,k1)',numSens,1) - sensors;
        Yk1k(:,k1) = sqrt(sum(xydiff.^2,2));
    end
    
    W = [lambda 0.5*ones(1,2*L)]/(L+lambda);
    Wmat = repmat(W,length(meas),1);
    
    ykminus = sum(Wmat.*Yk1k,2);
    
    Wmat(:,1) = Wmat(:,1) + (1-alpha^2+mybeta);
    
    ykminusmat = repmat(ykminus,1,num_sigmaPts);
    Pykyk = Wmat .* (Yk1k - ykminusmat)*(Yk1k - ykminusmat)';
    Pxkyk = (sigmaPts - xkminusmat)*(Wmat .*(Yk1k - ykminusmat))';
    
    Rk = Parameters.VarianceOfRangeMeasurement * eye(numSens);
    
    Pykyk = Pykyk + Rk;
    
    Kk = Pxkyk*inv(Pykyk);
    
    nuk = meas - ykminus;
    
    Xk1 = xkminus + Kk*nuk;
    Pk1 = Pkminus - Kk*Pykyk*Kk';
    
    likeli = (1/sqrt(det(2*pi*Pykyk)))*exp(-0.5*nuk'*inv(Rk)*nuk);
end


