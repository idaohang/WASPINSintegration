


% Plots the filter etimated sensor biases vs. the original data (call the
% main function first)

global  Pk timestamps_IMU nS 

scrsz = get(0,'ScreenSize'); % Full screen
figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)]);

subplot(3, 1, 1), plot(timestamps_IMU(1:nS),Pk(1:3,1:nS)), grid,
xlabel('time in secs'),
ylabel('position error m'), legend('x','y','z');

subplot(3, 1, 2), plot(timestamps_IMU(1:nS),Pk(4:6,1:nS)), grid,
xlabel('time in secs'),
ylabel('velocity error ms'), legend('x','y','z');

subplot(3, 1, 3), plot(timestamps_IMU(1:nS),rad2deg(Pk(7:9,1:nS))), grid,
xlabel('time in secs'),
ylabel('orientation error deg'), legend('R','P','Y');


figure('OuterPosition',[1+scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);

subplot(2, 1, 1), plot(timestamps_IMU(1:nS),Pk(10:12,1:nS)), grid,
xlabel('time in secs'),
ylabel('gyro bias error deg/s'), legend('R','P','Y');

subplot(2, 1, 2), plot(timestamps_IMU(1:nS),Pk(13:15,1:nS)), grid,
xlabel('time in secs'),
ylabel('acc bias error ms2'), legend('x','y','z');