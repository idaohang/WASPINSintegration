
% Plots the std dev of orientation and acc errors

global Pk timestamps_WASP;

scrsz = get(0,'ScreenSize'); % Full screen
figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)]);

subplot(3, 1, 1), plot(timestamps_WASP, rad2deg(Pk(7:9,2:end))), grid, 
legend('R','P','Y'), 
ylabel('Std dev of orientation error deg'),
xlabel('time in secs');

subplot(3, 1, 2), plot(timestamps_WASP, Pk(10:12,2:end)), grid, 
legend('R','P','Y'), 
ylabel('Std dev of gyro bias deg/s'),
xlabel('time in secs');

subplot(3, 1, 3), plot(timestamps_WASP, Pk(13:15,2:end)), grid, 
legend('x','y','z'), 
ylabel('Std dev of accelerometer bias mg'),
xlabel('time in secs');