function [Xk1, Pk1, innov, Sk_diag] = KF(xkm, Pkm, zk, Fk, Qk, Hk, Rk)
% Kalman Filter algorithm for a linear Gaussian state space model

% outputs of the function are:
% Xk1 - updated state
% Pk1 – updated state covariance

% inputs to the function are:
% xkm – prior state
% Pkm – prior state covariance
% zk – measurement
% Fk – state transition matrix
% Qk – process noise covariance
% Hk - measurement matrix
% Rk – measurement noise covariance

% Predicted state 
xk = Fk * xkm;
% predicted state covariance
Pk = Fk * Pkm * Fk' + Qk;
% predicted measurement
nuxk = Hk * xk;
% innovation
innov = zk - nuxk;
% predicted measurement covariance
Sk = Hk * Pk * Hk' + Rk;
% Kalman gain
Kk = Pk * Hk' / Sk;
% update estimate with measurement
Xk1 = xk + Kk*(innov);
% compute error covariance for updated estimate
Pk1 = (eye(size(Kk*Hk)) - Kk*Hk) * Pk;

Sk_diag = sqrt(diag(Sk));

end

