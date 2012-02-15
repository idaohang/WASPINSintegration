


% Plots the position, velocity, and orientation corrections estimated by
% the filter (call the main function first)

global dr dv epsilon timestamps_WASP nS_filter 

scrsz = get(0,'ScreenSize'); % Full screen
figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)]);

subplot(3, 1, 1), plot(timestamps_WASP(1:nS_filter), dr(:,1:nS_filter), '.'), grid
xlabel('time in secs'),
ylabel('Est. Position error m') , legend('x','y','z');

subplot(3, 1, 2), plot(timestamps_WASP(1:nS_filter), dv(:,1:nS_filter),'.'), grid
xlabel('time in secs'),
ylabel('Est. Velocity error m/s') , legend('x','y','z');

subplot(3, 1, 3), plot(timestamps_WASP(1:nS_filter), rad2deg(epsilon(:,1:nS_filter)),'.'), grid
xlabel('time in secs'),
ylabel('Est. Orientation error epsilon deg') , legend('R','P','Y');


%figure('OuterPosition',[1+scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);







