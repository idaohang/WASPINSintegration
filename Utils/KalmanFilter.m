% An implementation of the Kalman Filter. Assumes nearly constant velocity
% model for state transition and position measurements.
%
% Input:
%   Xk         - Prior state.
%   Pk         - Prior state covariance.
%   zk         - Position measurement.
%   Parameters - A structure consisting of various parameters.
%
% Output:
%   Xk1        - Updated state.
%   Pk1        - Updated state covariance

function [Xk1,Pk1] = KalmanFilter(Xk, Pk, zk, Parameters)

[Fk,Qk] = GetMotionModel(Parameters.SpaceDimension, 'StateTransitionMatrix_CV', ...
    Parameters.SamplingTime, 'ProcessNoiseCovarianceMatrix_CV', Parameters.SamplingTime, Parameters.ProcessNoiseIntensity);

Hk = [1 0];
if (Parameters.SpaceDimension == 2)
    Hk = blkdiag(Hk,Hk);
elseif (Parameters.SpaceDimension == 3)
    Hk = blkdiag(Hk,Hk,Hk);
end

if (Parameters.SpaceDimension == 2)
    Rk = blkdiag(Parameters.VarianceOfPositionEstimates, Parameters.VarianceOfPositionEstimates);
elseif (Parameters.SpaceDimension == 3)
    Rk = blkdiag(Parameters.VarianceOfPositionEstimates, Parameters.VarianceOfPositionEstimates, Parameters.VarianceOfPositionEstimates);
end

%% uncomment the following lines to use CRLB as measurement covariance
%     if ~isempty(X1)
%         Rk = Calc_CRLB(X1, sensors, DIM) / .01;
%     else
%         Rk = [];
%     end

xkm = Xk;
Pkm = Pk;

xk = Fk * xkm;
Pk = Fk * Pkm * Fk' + Qk;

if isempty(zk)
    % no measurement to update
    Xk1 = xk;
    Pk1 = Pk;
else
    % have measurement - run normally
    % predicted measurement
    nuk = Hk * xk;
    % predicted measurement covariance
    Sk = Hk * Pk * Hk' + Rk;
    % compute Kalman gain
    Kk = Pk * Hk' * inv(Sk);
    
    % update estimate with measurement
    Xk1 = xk + Kk*(zk - nuk);
    % compute error covariance for updated estimate
    Pk1 = (eye(size(Kk*Hk)) - Kk*Hk) * Pk * (eye(size(Kk*Hk)) - Kk*Hk) + Kk * Rk * Kk';
end

