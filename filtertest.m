



% Testing the errorfilter.m function
clear all;

Cs2b = eye(3);
Cb2n = (RPY2DCM([deg2rad(0);deg2rad(0);deg2rad(0)]))' % Sitting upright in n-frame at 0, -90, 0
error = (RPY2DCM([deg2rad(-2);deg2rad(3);deg2rad(-2)]))'
dt = 1; % one second for easier to work with

dr_wasp = [0; 0; 0] ; % some x error in nav frame, is -z for sensor

%mn_cal = [209; 113; 450];
mn_cal = [1; 0; 0];
gn_cal = [0; 0; -1];

%ms = Cb2n'*[0.7; 0.7; 0.7] % less x, more y, = -yaw bias
ms = Cb2n' * mn_cal
fs = Cb2n' * gn_cal * 1000 % gravity only

Cb2n_error = error*Cb2n;


T=5;
for t=1:T
    [ ba, bg, dv, dr, epsilon, innov, Sk_diag, Pk_diag ]...
        = errorfilter( Cs2b, Cb2n_error, fs, ms, dt, dr_wasp, mn_cal, gn_cal );
    Cb2n_adjusted = (eye(3) + mat_cross(epsilon)) * Cb2n_error
    epsilon = rad2deg(epsilon)
    Cb2n_error = Cb2n_adjusted;
end;
% Cb2n = (eye(3) + mat_cross(epsilon )) * Cb2n;
% 
% [ ba, bg, dv, dr, epsilon, innov, Sk_diag, Pk_diag ]...
%     = errorfilter( Cs2b, Cb2n, fs, ms, dt, dr_wasp, mn_cal, gn_cal );
% 
% epsilon = rad2deg(epsilon)
% bg = rad2deg(bg)
