


% Plots the filter etimated sensor biases vs. the original data (call the
% main function first)

global  fs ba bg timestamps_IMU timestamps_WASP nS nS_filter actual_ba actual_bg fsc ws wsc

scrsz = get(0,'ScreenSize'); % Full screen
figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)]);

subplot(4, 1, 1), plot(timestamps_IMU(1:nS),fs(:,1:nS)), grid
xlabel('time in secs'),
ylabel('Raw Accel mg'), legend('x','y','z');

subplot(4, 1, 2), plot(timestamps_WASP(1:nS_filter),ba(:,1:nS_filter)), grid
xlabel('time in secs'),
ylabel('KF Accel bias mg'), legend('x','y','z');

% Correct if stationary only!
actual_ba_plot = MovAvg2(820, actual_ba')';
subplot(4, 1, 3), plot(timestamps_IMU(1:nS),actual_ba_plot(:,1:nS)), grid
xlabel('time in secs'),
ylabel('Actual Accel bias (when static) mg'), legend('x','y','z');

subplot(4, 1, 4), plot(timestamps_IMU(1:nS),fsc(:,1:nS)), grid
xlabel('time in secs'),
ylabel('Corrected Accel mg'), legend('x','y','z');



figure('OuterPosition',[1+scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);

subplot(4, 1, 1), plot(timestamps_IMU(1:nS),ws(:,1:nS)), grid
xlabel('time in secs'),
ylabel('Raw gyro deg/s'), legend('R','P','Y');

subplot(4, 1, 2), plot(timestamps_WASP(1:nS_filter),bg(:,1:nS_filter)), grid
xlabel('time in secs'),
ylabel('Gyro bias deg/s') , legend('R','P','Y');

% Correct if stationary only!
actual_bg_plot = MovAvg2(1600, actual_bg')';
subplot(4, 1, 3), plot(timestamps_IMU(1:nS),actual_bg_plot(:,1:nS)), grid
xlabel('time in secs'),
ylabel('Actual Gyro bias (when static) deg/s'), legend('x','y','z');


subplot(4, 1, 4), plot(timestamps_IMU(1:nS),wsc(:,1:nS)), grid
xlabel('time in secs'),
ylabel('Corrected gyro deg/s'), legend('x','y','z');





