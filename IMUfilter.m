function [ bg, epsilon, innov, Sk_diag, Pkm]...
    = IMUfilter( Cs2b, Cb2n, fs, ms, dt, mn_cal, gn_cal, fs_var, Pkm )
%errorfilter Implements an error state filter
%   outputs of the function are:
%   bg - estimated gyro bias in the s-frame (in deg/s)
%   epsilon - estimated orientation error in Euler angles around the
%       x, y, z axes respectively (roll, pitch, yaw) in radians

%   inputs to the function are:
%   Cs2b - constant Direction Cosine Matrix (DCM) from the s-frame to the b-frame
%   Cb2n - estimated DCM from the b-frame to the n-frame
%   fs - bias compensated accelerometer measurement in the s-frame (in mg)
%   ms - magnetometer measurement in the s-frame (in mgauss)
%   dt - sampling interval in seconds
%   mn_cal - calibrated environment magnetic field in the n-frame\
%   gn_cal - calibrated gravity in the n-frame
%   fs_var - the variance of the accelerometer data (in mg squared)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FILTER TUNING PARAMETERS

%   Initial state covariance
Pinitial = diag([0.01; 0.01; 0.01; 1; 1; 1]);

%   Sensors noise covariances
Ng = 0.1;   % Gyro sensor noise cov (in deg/s squared) 1 % increase to make bias estimate more stable
Nm = 500;    % Magnetometer sensor noise covariance 50

Ug = 0.005;   % Gyro bias drift noise cov 0.01 % decrease to make bias estimate stable

%   Adaptive covariance weighting scaling factor
alpha_f = 5;          % Accel observation 20 100
alpha_m = 30;          % Magnetometer observation 30 100


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   MATRIX CONSTRUCTION FOR KALMAN FILTER
%   Refer to WASP / INS integration paper for further details

%   Temporary variables
Cs2n = Cb2n*Cs2b;           % estimated DCM s-frame to n-frame
Cn2s = Cs2n';               % estimated DCM n-frame to s-frame

gn_mag = norm(gn_cal);     
gn_dir = gn_cal/gn_mag;
fs_mag = norm(fs);          % magnitude of accel in s-frame (in mg)
fs_dir = fs/fs_mag;         % direction of accel in s-frame (unit length)       
ms_mag = norm(ms);          % magnitude of magnetic field in s-frame (in mgauss)
ms_dir = ms/ms_mag;         % direction of magnetic in s-frame (unit length)
mn_cal_mag = norm(mn_cal);        % magnitude of magnetic calibration in n-frame (in mgauss)
mn_cal_dir = mn_cal/mn_cal_mag;   % direction of magnetic calibration in n-frame (unit length)
  

%   Prior state - always reset to zero, since we feedback errors every time
%   orientation error, gyro bias error, magnetometer bias
xkm = zeros(6,1); 

%   State transition matrix: orientation error, gyro bias, magnetometer bias 
Ft = [  zeros(3,3),     deg2rad(-Cs2n);
        zeros(3,3),     zeros(3,3)];

Fk = expm(Ft*dt);    % Can speed this up with a lower order approximation later on

%   Prior state covariance
%persistent Pkm;
if isempty(Pkm)
    Pkm = Pinitial;
end;

%   Process noise covariance
Qt = diag([Ng*ones(3,1); Ug*ones(3,1)]);

G = [  deg2rad(-Cs2n), zeros(3,3);
       zeros(3,3),     eye(3)];

Qk = 1/2*(Fk*G*Qt*G' + G*Qt*G'*Fk')*dt;


%   Measurement noise covariance (adaptive) refer "Adaptive Sensor Data
%   Fusion in Motion Capture"
ms_mag_error = max([ ms_mag/mn_cal_mag; mn_cal_mag/ms_mag]); % from 1 to Inf
fs_mag_error = max([ fs_mag/1000; 1000/fs_mag]);   % fs magnitude should be 1000mg

ms_dir_error = norm( Cs2n*ms_dir - mn_cal_dir ); % from 0 to Inf
fs_dir_error = norm( Cs2n*fs_dir - gn_dir );

ms_confidence = exp(alpha_m*ms_mag_error*ms_dir_error); % exponent is from 0 to Inf
fs_confidence = exp(alpha_f*fs_mag_error*fs_dir_error);

%   Measurement
ms_error = Cs2n*ms - mn_cal;
ms_error(3)=0;
zk = [ fs - Cn2s*gn_cal;      Cn2s*ms_error]; 

%zk = [ fs - Cn2s*gn_cal;      ms - Cn2s*mn_cal];    

%   Measurement matrix
Hk = [ Cn2s*mat_cross(gn_cal),      zeros(3,3); 
       Cn2s*mat_cross(mn_cal),      zeros(3,3)]; 

%   Measurement covariance
Rk = diag([min(fs_var*ones(3,1)*fs_confidence, 5000^2);  min(Nm*ones(3,1)*ms_confidence, 100^2) ]);
%   (need to cap variance for precision of matrix)

%   Run the Kalman Filter
[xkm, Pkm, innov, Sk_diag] = KF(xkm, Pkm, zk, Fk, Qk, Hk, Rk);
epsilon = xkm(1:3,:);
bg = xkm(4:6,:);

end




