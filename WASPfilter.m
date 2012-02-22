function [ ba, bg, dv, dr, epsilon, innov, Sk_diag, Pk_diag ]...
    = WASPfilter( Cs2b, Cb2n, fs_mean, dt, dr_wasp )
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
%   fs_mean - bias compensated accelerometer measurement in the s-frame (in mg)
%   dt - sampling interval in seconds
%   dr_wasp - observed INS position error based on WASP (in m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FILTER TUNING PARAMETERS

%   Initial state covariance
Pinitial = diag([2; 2; 2; 2; 2; 2; 0.01; 0.01; 0.01; 1; 1; 1; 100; 100; 100]);

%   WASP covariance
wasp_err = 0.10^2; % (in m) 0.10
vert_err = 0.5^2; % Assume vertical movement is an error, but not much certainty

%   Sensors noise covariances
Na = 100; % Accel sensor noise cov (in mg squared) 100 % increase to put less weight on observation
Ng = 1;   % Gyro sensor noise cov (in deg/s squared) 1 % increase to make bias estimate more stable

Ua = 1;  % Accel bias drift noise cov 30
Ug = 0;   % Gyro bias drift noise cov 0.01 % decreas to make bias estimate stable


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   MATRIX CONSTRUCTION FOR KALMAN FILTER
%   Refer to WASP / INS integration paper for further details

%   Temporary variables
Cs2n = Cb2n*Cs2b;           % estimated DCM s-frame to n-frame
fn_cross = mg2ms(mat_cross(Cs2n*fs_mean));


%   Prior state - always reset position, velocity, orientation error to zero, since we feedback errors every time
xkm = zeros(15,1); 

%   State transition matrix     
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
zk = [ dr_wasp ];  

%   Measurement matrix
Hk = [ eye(3),     zeros(3,3),     zeros(3,3),                  zeros(3,3),     zeros(3,3)];   

%   Measurement covariance
Rk = diag([wasp_err*ones(2,1); vert_err]);
 

%   Run the Kalman Filter
[xkm, Pkm, innov, Sk_diag] = KF(xkm, Pkm, zk, Fk, Qk, Hk, Rk);
dr = xkm(1:3,:);
dv = xkm(4:6,:);
epsilon = xkm(7:9,:);
bg = xkm(10:12,:);
ba = xkm(13:15,:);

Pk_diag = sqrt(diag(Pkm));

end




