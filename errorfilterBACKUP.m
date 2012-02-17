function [ ba, bg, dv, dr, epsilon, innov, Sk_diag, Pk_diag ]...
    = errorfilter( Cs2b, Cb2n, fs, ms, dt, dr_wasp, mn_cal, gn_cal )
%errorfilter Implements an error state filter
%   outputs of the function are:
%   ba - estimated accelerometer bias in the s-frame (in mg)
%   bg - estimated gyro bias in the s-frame (in deg/s)
%   dv - estimated velocity error in the n-frame (in m/s)
%   dr - estimated position error in the n-frame (in m)
%   epsilon - estimated orientation error in Euler angles around the
%       x, y, z axes respectively (roll, pitch, yaw) in radians

%   inputs to the function are:
%   Cs2b - constant Direction Cosine Matrix (DCM) from the s-frame to the b-frame
%   Cb2n - estimated DCM from the b-frame to the n-frame
%   fs - bias compensated accelerometer measurement in the s-frame (in mg)
%   ms - magnetometer measurement in the s-frame (in mgauss)
%   dt - sampling interval in seconds
%   dr_wasp - observed INS position error based on WASP (in m)
%   wasp_err - a measure of the error in the observed WASP location
%   mn_cal - calibrated environment magnetic field in the n-frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FILTER TUNING PARAMETERS

%   Initial state covariance
Pinitial = diag([2; 2; 2; 2; 2; 2; 0.1; 0.1; 0.1; 0.5; 0.5; 0.5; 100; 100; 100]);

%   WASP covariance
wasp_err = 0.10^2; % (in m) 0.10
vert_err = 0.5^2; % Assume vertical movement is an error, but not much certainty

%   Sensors noise covariances
Na = 100;    % Accel sensor noise cov (in mg squared) 50 % increase to put less weight on observation
Ng = 1;   % Gyro sensor noise cov (in deg/s squared) 2 % increase to make bias estimate more stable
Nm = 25;    % Magnetometer sensor noise covariance

Ua = 1;  % Accel bias drift noise cov 30
Ug = 0;   % Gyro bias drift noise cov 0.01 % decreas to make bias estimate stable

%Na_obs =  0.0005;  % Noise in each component of the accel and mag direction vector 0.0005 0.00005
%Nm_obs =  0.00005;  % Noise in each component of the accel and mag direction vector 0.00005
% Experience on a stationary node suggests 0.0001, therefore use scaling
% factor to increase, not this
 
%Ta = 10;     % Accel correlation time in seconds (close 0 = random walk) 10
%Tg = 3;  % Gyro correlation time in intervals (close to 0 = random walk)

%   Adaptive covariance weighting scaling factor
alpha_f = 50;            % Accel observation 20 100
alpha_m = 0;          % Magnetometer observation 30 100



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   MATRIX CONSTRUCTION FOR KALMAN FILTER
%   Refer to WASP / INS integration paper for further details

%   Temporary variables
Cs2n = Cb2n*Cs2b;           % estimated DCM s-frame to n-frame
Cn2s = Cs2n';               % estimated DCM n-frame to s-frame

gn_mag = norm(gn_cal);     
gn_dir = gn_cal/gn_mag;

fn = (Cs2n*fs - gn_dir*1000);   % accel in n-frame excluding gravity 
%fn_cross = mg2ms(mat_cross(fn)); % matrix form and units conversion to m/s2

fn_cross = mg2ms(mat_cross(Cs2n*fs));

fs_mag = norm(fs);          % magnitude of accel in s-frame (in mg)
fs_dir = fs/fs_mag;         % direction of accel in s-frame (unit length)       
ms_mag = norm(ms);          % magnitude of magnetic field in s-frame (in mgauss)
ms_dir = ms/ms_mag;         % direction of magnetic in s-frame (unit length)
mn_cal_mag = norm(mn_cal);        % magnitude of magnetic calibration in n-frame (in mgauss)
mn_cal_dir = mn_cal/mn_cal_mag;   % direction of magnetic calibration in n-frame (unit length)
  

%   Prior state - always reset position, velocity, orientation error to zero, since we feedback errors every time
xkm = zeros(15,1); 
%   State transition matrix - note that the Time constants can be put to
%   diag((-1/Tg)*ones(3,1))
    
Ft = [  zeros(3,3), eye(3),     zeros(3,3),     zeros(3,3),        zeros(3,3);
        zeros(3,3), zeros(3,3), fn_cross,       zeros(3,3),        mg2ms(Cs2n);
        zeros(3,3), zeros(3,3), zeros(3,3),     deg2rad(-Cs2n),    zeros(3,3);
        zeros(3,3), zeros(3,3), zeros(3,3),     zeros(3,3),        zeros(3,3);
        zeros(3,3), zeros(3,3), zeros(3,3),     zeros(3,3),        zeros(3,3)];


Fk = expm(Ft*dt);    % Can speed this up with a lower order approximation later on

%   Prior state covariance
persistent Pkm;
if isempty(Pkm)
    Pkm = Pinitial;
end;

%   Process noise covariance
Qt = diag([Na*ones(3,1); Ng*ones(3,1); Ug*ones(3,1); Ua*ones(3,1)]);

G = [  zeros(3,3),          zeros(3,3),     zeros(3,3),      zeros(3,3);
       mg2ms(Cs2n),         zeros(3,3),     zeros(3,3),      zeros(3,3);
       zeros(3,3),          deg2rad(-Cs2n), zeros(3,3),      zeros(3,3);
       zeros(3,3),          zeros(3,3),     eye(3),          zeros(3,3);
       zeros(3,3),          zeros(3,3),     zeros(3,3),      eye(3)];

Qk = 1/2*(Fk*G*Qt*G' + G*Qt*G'*Fk')*dt;

%   Measurement
zk = [ dr_wasp;     fs - Cn2s*gn_cal;      ms - Cn2s*mn_cal];   % changed to the n-frame for magnetometer

% mag_error =  Cs2n*ms_dir - mn_cal_dir;
% mag_error(3) = 0;    % can use mag only for yaw error
% zk = [ dr_wasp;           fs_dir - Cn2s*gn_dir;          mag_error];

% mag_error =  Cs2n*ms - mn_cal;
% mag_error(3) = 0;    % can use mag only for yaw error
% zk = [ dr_wasp;           fs - Cn2s*gn_cal;          mag_error];

%   Measurement matrix
% Hk = [ eye(3),     zeros(3,3),     zeros(3,3),                  zeros(3,3),     zeros(3,3);
%        zeros(3,3), zeros(3,3),     Cn2s*mat_cross(gn_dir),      zeros(3,3),     zeros(3,3); 
%        zeros(3,3), zeros(3,3),     mat_cross(mn_cal_dir),       zeros(3,3),     zeros(3,3)];
  
Hk = [ eye(3),     zeros(3,3),     zeros(3,3),                  zeros(3,3),     zeros(3,3);
       zeros(3,3), zeros(3,3),     Cn2s*mat_cross(gn_cal),      zeros(3,3),     eye(3); 
       zeros(3,3), zeros(3,3),     Cn2s*mat_cross(mn_cal),      zeros(3,3),     zeros(3,3)];
   
%Hk = [  eye(3),     zeros(3,3),     zeros(3,3),                zeros(3,3),     zeros(3,3)];

%   Measurement noise covariance (adaptive) refer "Adaptive Sensor Data
%   Fusion in Motion Capture"
ms_mag_error = max([ ms_mag/mn_cal_mag; mn_cal_mag/ms_mag]); % from 1 to Inf
fs_mag_error = max([ fs_mag/1000; 1000/fs_mag]);   % fs magnitude should be 1000mg

ms_dir_error = norm( Cs2n*ms_dir - mn_cal_dir ); % from 0 to Inf
fs_dir_error = norm( Cs2n*fs_dir - gn_dir );

ms_confidence = exp(alpha_m*ms_mag_error*ms_dir_error); % exponent is from 0 to Inf
fs_confidence = exp(alpha_f*fs_mag_error*fs_dir_error);

% Rk = diag([wasp_err*ones(2,1); vert_err; min(fs_confidence*Na_obs,10^2)*ones(3,1);  min(ms_confidence*Nm_obs,10^2)*ones(3,1)]);

Rk = diag([wasp_err*ones(2,1); vert_err; min(fs_confidence*Na,1000^2)*ones(3,1);  min(ms_confidence*Nm,100^2)*ones(3,1)]);

% Rk = diag([wasp_err*ones(2,1); vert_err;                     min(ms_confidence*Nm_obs,10^2)*ones(3,1)]);

% if (ms_confidence*Na_obs > 10^2)
%     zk(4:6) = [0; 0; 0];   
% end;

%   Run the Kalman Filter
[xkm, Pkm, innov, Sk_diag] = KF(xkm, Pkm, zk, Fk, Qk, Hk, Rk);
dr = xkm(1:3,:);
dv = xkm(4:6,:);
epsilon = xkm(7:9,:);
bg = xkm(10:12,:);
ba = xkm(13:15,:);

Pk_diag = sqrt(diag(Pkm));

% For modifying observations
% innov = [innov(1:3); zeros(3,1); innov(4:6)];
% Sk_diag = [Sk_diag(1:3); zeros(3,1); Sk_diag(4:6)];

end




