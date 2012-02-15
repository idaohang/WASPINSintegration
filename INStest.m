



% Testing the errorfilter.m function
clear all;

Cb2n = (RPY2DCM([deg2rad(0);deg2rad(-90);deg2rad(0)]))'; % Sitting upright in n-frame
wb=deg2rad([45;0;0]);
dt=0.1;

Cb2n'

%Cb2n =  Cb2n * expm(mat_cross(wb*dt));
for i=1:10
    Cb2n =  (Cb2n' * (eye(3) + mat_cross(wb*dt)))';
end;

Cb2n'
x = Cb2n'*[1;0;0]
y = Cb2n'*[0;1;0]