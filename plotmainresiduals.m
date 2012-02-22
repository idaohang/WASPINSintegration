

% Plots the filter observation residuals and dynamically scaled error bands
% (call the main function first)

global timestamps_IMU nS Sk innov


% Plot position residuals
scrsz = get(0,'ScreenSize'); % Full screen
figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)]);

subplot(6, 1, 1), plot( timestamps_IMU(1:nS), innov(1,1:nS), 'b-' ), grid
ylabel('X pos error innov m'),
xlabel('time in secs'),
%axis([0; timestamps_IMU(nS); xmin; xmax]),
hold on,
plot(timestamps_IMU(1:nS), Sk(1,1:nS), 'r--'), grid
plot(timestamps_IMU(1:nS), -Sk(1,1:nS), 'r--'), grid
hold off;

subplot(6, 1, 2), plot( timestamps_IMU(1:nS), innov(2,1:nS), 'b-' ), grid
ylabel('Y pos error innov m'),
xlabel('time in secs'),
%axis([0; timestamps_IMU(nS); xmin; xmax]),
hold on,
plot(timestamps_IMU(1:nS), Sk(2,1:nS), 'r--'), grid
plot(timestamps_IMU(1:nS), -Sk(2,1:nS), 'r--'), grid
hold off;

subplot(6, 1, 3), plot( timestamps_IMU(1:nS), innov(3,1:nS), 'b-' ), grid
ylabel('Z pos error innov m'),
xlabel('time in secs'),
%axis([0; timestamps_IMU(nS); xmin; xmax]),
hold on,
plot(timestamps_IMU(1:nS), Sk(3,1:nS), 'r--'), grid
plot(timestamps_IMU(1:nS), -Sk(3,1:nS), 'r--'), grid
hold off;

% Change scaling

xmin = -0.5;
xmax = 0.5;

% Plot gravity error residuals

subplot(6, 1, 4), plot( timestamps_IMU(1:nS), innov(4,1:nS), 'b-' ), grid
ylabel('X Gravity error innov'),
xlabel('time in secs'),
%axis([0; timestamps_IMU(nS); xmin; xmax]),
hold on,
plot(timestamps_IMU(1:nS), Sk(4,1:nS), 'r--'), grid
plot(timestamps_IMU(1:nS), -Sk(4,1:nS), 'r--'), grid
hold off;

subplot(6, 1, 5), plot( timestamps_IMU(1:nS), innov(5,1:nS), 'b-' ), grid
ylabel('Y Gravity error innov'),
xlabel('time in secs'),
%axis([0; timestamps_IMU(nS); xmin; xmax]),
hold on,
plot(timestamps_IMU(1:nS), Sk(5,1:nS), 'r--'), grid
plot(timestamps_IMU(1:nS), -Sk(5,1:nS), 'r--'), grid
hold off;

subplot(6, 1, 6), plot( timestamps_IMU(1:nS), innov(6,1:nS), 'b-' ), grid
ylabel('Z Gravity error innov'),
xlabel('time in secs'),
%axis([0; timestamps_IMU(nS); xmin; xmax]),
hold on,
plot(timestamps_IMU(1:nS), Sk(6,1:nS), 'r--'), grid
plot(timestamps_IMU(1:nS), -Sk(6,1:nS), 'r--'), grid
hold off;

% Plot magnetic error residuals
figure('OuterPosition',[1+scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);

subplot(3, 1, 1), plot( timestamps_IMU(1:nS), innov(7,1:nS), 'b-' ), grid
ylabel('X Magnetic error innov'),
xlabel('time in secs'),
%axis([0; timestamps_IMU(nS); xmin; xmax]),
hold on,
plot(timestamps_IMU(1:nS), Sk(7,1:nS), 'r--'), grid
plot(timestamps_IMU(1:nS), -Sk(7,1:nS), 'r--'), grid
hold off;

subplot(3, 1, 2), plot( timestamps_IMU(1:nS), innov(8,1:nS), 'b-' ), grid
ylabel('Y Magnetic error innov'),
xlabel('time in secs'),
%axis([0; timestamps_IMU(nS); xmin; xmax]),
hold on,
plot(timestamps_IMU(1:nS), Sk(8,1:nS), 'r--'), grid
plot(timestamps_IMU(1:nS), -Sk(8,1:nS), 'r--'), grid
hold off;

subplot(3, 1, 3), plot( timestamps_IMU(1:nS), innov(9,1:nS), 'b-' ), grid
ylabel('Z Magnetic error innov'),
xlabel('time in secs'),
%axis([0; timestamps_IMU(nS); xmin; xmax]),
hold on,
plot(timestamps_IMU(1:nS), Sk(9,1:nS), 'r--'), grid
plot(timestamps_IMU(1:nS), -Sk(9,1:nS), 'r--'), grid
hold off;




