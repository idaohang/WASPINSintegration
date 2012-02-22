function [ a, v, r, Cb2n ] = INS( v, r, Cb2n, dv, dr, epsilon, wb, fb, dt, g)
%INS implements an Inertial Navigation System update
%   a is the body acceleration in the n-frame (in m/s^2)
%   v is the body velocity in the n-frame (in m/s)
%   r is the body position in the n-frame (in m)
%   Cb2n is the DCM from the b-frame to the n-frame
%   dv is the estimated velocity error in the n-frame (in m/s)
%   dr is the estimated position error in the n-frame (in m)
%   epsilon - estimated orientation error in Euler angles around the
%       x, y, z axes respectively (roll, pitch, yaw) in radians
%   wb is the bias corrected gyro in the b-frame (in deg/s)
%   fb is the bias corrected accelerometer in the b-frame (in mg)
%   dt is the time increment 

wb = deg2rad(wb);

% Acceleration update (including the rotation correction)
a = mg2ms(Cb2n * (fb  + 0.5 * cross(wb*dt, fb)) - g);

% Velocity update with error correction
v_last = v;
v = v  + a*dt - dv;

% Position update by trapezoidal integration with error correction
r = r + 0.5*(v + v_last)*dt - dr;

% Orientation update with error correction
Cb2n = (eye(3) + mat_cross(epsilon)) * Cb2n;
Cb2n = (Cb2n * (eye(3) + mat_cross(wb*dt)));

% Now re-orthogonalise and re-normalize the orientation matrix (at a lower frequency)
persistent count;
if isempty(count)
    count = 0;
end;
count = count + 1;

if count == 50 % Every 50 iterations
    count = 0;
    [U S V] = svd(Cb2n);
    Cb2n = U*V'; % Drop the diagonal
end

end

