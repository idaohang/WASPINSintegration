

% Plots the filter observation residuals and dynamically scaled error bands
% (call the main function first)

global timestamps_WASP nS_filter Sk innov


% Plot position residuals
scrsz = get(0,'ScreenSize'); % Full screen
figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)]);

subplot(6, 1, 1), plot( timestamps_WASP(1:nS_filter), innov(1,1:nS_filter), 'b-' ), grid
ylabel('X pos error innov m'),
xlabel('time in secs'),
%axis([0; timestamps_WASP(nS_filter); xmin; xmax]),
hold on,
plot(timestamps_WASP(1:nS_filter), Sk(1,1:nS_filter), 'r--'), grid
plot(timestamps_WASP(1:nS_filter), -Sk(1,1:nS_filter), 'r--'), grid
hold off;

subplot(6, 1, 2), plot( timestamps_WASP(1:nS_filter), innov(2,1:nS_filter), 'b-' ), grid
ylabel('Y pos error innov m'),
xlabel('time in secs'),
%axis([0; timestamps_WASP(nS_filter); xmin; xmax]),
hold on,
plot(timestamps_WASP(1:nS_filter), Sk(2,1:nS_filter), 'r--'), grid
plot(timestamps_WASP(1:nS_filter), -Sk(2,1:nS_filter), 'r--'), grid
hold off;

subplot(6, 1, 3), plot( timestamps_WASP(1:nS_filter), innov(3,1:nS_filter), 'b-' ), grid
ylabel('Z pos error innov m'),
xlabel('time in secs'),
%axis([0; timestamps_WASP(nS_filter); xmin; xmax]),
hold on,
plot(timestamps_WASP(1:nS_filter), Sk(3,1:nS_filter), 'r--'), grid
plot(timestamps_WASP(1:nS_filter), -Sk(3,1:nS_filter), 'r--'), grid
hold off;

% Change scaling

xmin = -0.5;
xmax = 0.5;

% Plot gravity error residuals

subplot(6, 1, 4), plot( timestamps_WASP(1:nS_filter), innov(4,1:nS_filter), 'b-' ), grid
ylabel('X Gravity error innov'),
xlabel('time in secs'),
%axis([0; timestamps_WASP(nS_filter); xmin; xmax]),
hold on,
plot(timestamps_WASP(1:nS_filter), Sk(4,1:nS_filter), 'r--'), grid
plot(timestamps_WASP(1:nS_filter), -Sk(4,1:nS_filter), 'r--'), grid
hold off;

subplot(6, 1, 5), plot( timestamps_WASP(1:nS_filter), innov(5,1:nS_filter), 'b-' ), grid
ylabel('Y Gravity error innov'),
xlabel('time in secs'),
%axis([0; timestamps_WASP(nS_filter); xmin; xmax]),
hold on,
plot(timestamps_WASP(1:nS_filter), Sk(5,1:nS_filter), 'r--'), grid
plot(timestamps_WASP(1:nS_filter), -Sk(5,1:nS_filter), 'r--'), grid
hold off;

subplot(6, 1, 6), plot( timestamps_WASP(1:nS_filter), innov(6,1:nS_filter), 'b-' ), grid
ylabel('Z Gravity error innov'),
xlabel('time in secs'),
%axis([0; timestamps_WASP(nS_filter); xmin; xmax]),
hold on,
plot(timestamps_WASP(1:nS_filter), Sk(6,1:nS_filter), 'r--'), grid
plot(timestamps_WASP(1:nS_filter), -Sk(6,1:nS_filter), 'r--'), grid
hold off;

% Plot magnetic error residuals
figure('OuterPosition',[1+scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);

subplot(3, 1, 1), plot( timestamps_WASP(1:nS_filter), innov(7,1:nS_filter), 'b-' ), grid
ylabel('X Magnetic error innov'),
xlabel('time in secs'),
%axis([0; timestamps_WASP(nS_filter); xmin; xmax]),
hold on,
plot(timestamps_WASP(1:nS_filter), Sk(7,1:nS_filter), 'r--'), grid
plot(timestamps_WASP(1:nS_filter), -Sk(7,1:nS_filter), 'r--'), grid
hold off;

subplot(3, 1, 2), plot( timestamps_WASP(1:nS_filter), innov(8,1:nS_filter), 'b-' ), grid
ylabel('Y Magnetic error innov'),
xlabel('time in secs'),
%axis([0; timestamps_WASP(nS_filter); xmin; xmax]),
hold on,
plot(timestamps_WASP(1:nS_filter), Sk(8,1:nS_filter), 'r--'), grid
plot(timestamps_WASP(1:nS_filter), -Sk(8,1:nS_filter), 'r--'), grid
hold off;

subplot(3, 1, 3), plot( timestamps_WASP(1:nS_filter), innov(9,1:nS_filter), 'b-' ), grid
ylabel('Z Magnetic error innov'),
xlabel('time in secs'),
%axis([0; timestamps_WASP(nS_filter); xmin; xmax]),
hold on,
plot(timestamps_WASP(1:nS_filter), Sk(9,1:nS_filter), 'r--'), grid
plot(timestamps_WASP(1:nS_filter), -Sk(9,1:nS_filter), 'r--'), grid
hold off;




