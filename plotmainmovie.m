


% Plot the WASP movements on the map (large version)

scrsz = get(0,'ScreenSize'); % Full screen

map = loadmap;
set(map, 'OuterPosition', [1 1 scrsz(3) scrsz(4)]);
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


% Plot stuff
jump = 10;

% Prepare the new file.
%vidObj = VideoWriter('peaks.avi');
%open(vidObj);

front = [1; 0; 0];  % for trolley runs
%front = [0; 0; -1];  % for backpack runs
direct = Cn2b(:,:,1)*front

txt = text(xmin+1,ymin+1,'0.0 seconds','FontSize',20); % add timestamps
arrow = quiver(r(1,1),r(2,1),direct(1),direct(2),'wo','LineWidth',2); % show direction


legend('Anchor','Position','WASP observation');

fps = 5;

for j = 2:imu_hz/fps:nS
    
    plot(r(1,j-1:j),r(2,j-1:j),'g.') % Plot integrated INS position
    %scatter(WASPloc(1,j),WASPloc(2,j),'rx') % Plot WASP observations
    
    
    
    direct = Cn2b(:,:,j)*front;
    set(arrow, 'XData', r(1,j), 'YData', r(2,j), 'UData', direct(1), 'VData', direct(2));
    set(txt,'string',[num2str(timestamps(j)) ' seconds']);
    
    F(j) = getframe(gcf);
    
    % Write each frame to the file.
    %currFrame = getframe;
    %writeVideo(vidObj,currFrame);

    
end
movie(F,1)

% Close the file.
%close(vidObj);


%T = 9;
%text(WASPloc(1,1:T:nS_filter),WASPloc(2,1:T:nS_filter),num2str(timestamps_filter(1:T:nS_filter)')) % WASP timestamps
%J = jump*82;
%text(r(1,1:J:nS),r(2,1:J:nS),num2str(timestamps(1:J:nS)'), 'Color', [0; 0; 1]) % INS timestamps

hold off;


