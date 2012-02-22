


% Plots the filter etimated sensor biases vs. the original data (call the
% main function first)

global  ba bg bm sa timestamps_IMU timestamps_WASP nS nS_filter

scrsz = get(0,'ScreenSize'); % Full screen
figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)]);

num_plots = 2;
if size(bm,2) > 0 num_plots = num_plots + 1; end;
if size(sa,2) > 0 num_plots = num_plots + 1; end;

if size(ba,2) == nS+1
    subplot(num_plots, 1, 1), plot(timestamps_IMU(1:nS),ba(:,1:nS)), grid;
else
    subplot(num_plots, 1, 1), plot(timestamps_WASP(1:nS_filter),ba(:,1:nS_filter)), grid
end
hold on;
xlabel('time in secs'),
ylabel('KF Accel bias mg'), legend('x','y','z');
hold off;

if size(ba,2) == nS+1
    subplot(num_plots, 1, 2), plot(timestamps_IMU(1:nS),bg(:,1:nS)), grid;
else
    subplot(num_plots, 1, 2), plot(timestamps_WASP(1:nS_filter),bg(:,1:nS_filter)), grid
end
hold on;
xlabel('time in secs'),
ylabel('Gyro bias deg/s') , legend('R','P','Y');
hold off;

if size(bm,2) > 0
    subplot(num_plots, 1, 3), plot(timestamps_IMU(1:nS),bm(:,1:nS)), grid
    xlabel('time in secs'),
    ylabel('KF Mag bias mGauss'), legend('x','y','z');
end;

if size(sa,2) > 0
    subplot(num_plots, 1, num_plots), plot(timestamps_IMU(1:nS),sa(:,1:nS)), grid
    xlabel('time in secs'),
    ylabel('KF Accel scale factor'), legend('x','y','z');
end;



