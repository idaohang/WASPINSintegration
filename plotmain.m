


% Plot the WASP movements on the map (large version), as well as charts
% or the acceleration, velocity, position and rpy orientation (call the
% main function first)


global a v r rpy Cn2b WASPloc timestamps_IMU timestamps_WASP nS nS_filter front


scrsz = get(0,'ScreenSize'); % Full screen

map = loadmap;
set(map, 'OuterPosition', [1 1 scrsz(3)/2 scrsz(4)]);
hold on;

% Reset the axes on the map
xmin = round(min(WASPloc(1,:))-1);
xmax = round(max(WASPloc(1,:))+1);
ymin = round(min(WASPloc(2,:))-1);
ymax = round(max(WASPloc(2,:))+1);

if( (ymax-ymin)*1.6 < (xmax-xmin))
    extra = (xmax-xmin)/1.6-(ymax-ymin);
    prev = ylim;
    ymax = min(ymax + extra/2, prev(2));
    ymin = max(ymin - extra/2, prev(1));
else 
    extra = (ymax-ymin)*1.6-(xmax-xmin);
    prev = xlim;
    xmax = min(xmax + extra/2, prev(2));
    xmin = max(xmin - extra/2, prev(1));
end;

axis([xmin; xmax; ymin; ymax]),


jump = 10;
scatter(r(1,1:jump:nS),r(2,1:jump:nS),'g.') % Plot integrated INS position


% Plot orientation periodically - this is x axis in white
step = 410;
direct = zeros(3,size(1:step:nS,2));
j = 1;
for i=1:step:nS
    direct(:,j) = Cn2b(:,:,i)'*front;
    j = j+1;
end;
quiver(r(1,1:step:nS),r(2,1:step:nS),direct(1,:),direct(2,:),'mo','LineWidth',2, 'AutoScaleFactor', 0.2); % show direction

scatter(WASPloc(1,1:nS_filter),WASPloc(2,1:nS_filter),'rx') % Plot WASP observations



T = 9;
text(WASPloc(1,1:T:nS_filter),WASPloc(2,1:T:nS_filter),num2str(timestamps_WASP(1:T:nS_filter)')) % WASP timestamps
J = jump;%*40;
%text(r(1,1:J:nS),r(2,1:J:nS),num2str(timestamps_IMU(1:J:nS)'), 'Color', [0; 0; 1]) % INS timestamps

legend('WASP anchor node','IMU aided track','IMU orientation','WASP position');
%legend('WASP anchor node','IMU aided track','WASP position');
hold off;


% Plots the INS output after calling the main function first
figure('OuterPosition',[1+scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);

subplot(4, 1, 1), plot(timestamps_IMU(1:nS),a(:,1:nS)), grid
xlabel('time in secs'),
ylabel('Acc ms-2'), legend('x','y','z');
subplot(4, 1, 2), plot(timestamps_IMU(1:nS),v(:,1:nS)), grid
xlabel('time in secs'),
ylabel('Vel ms-1'), legend('x','y','z');
subplot(4, 1, 3), plot(timestamps_IMU(1:nS),r(:,1:nS)), grid
xlabel('time in secs'),
ylabel('Position m'), legend('x','y','z');
subplot(4, 1, 4), plot(timestamps_IMU(1:nS),rpy(:,1:nS)), grid
xlabel('time in secs'),
ylabel('Roll Pitch Yaw deg') , legend('R','P','Y');