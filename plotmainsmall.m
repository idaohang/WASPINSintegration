

% Plot the WASP movements on the map for just the first few seconds
% (call the main function first)


global r WASPloc timestamps_WASP

scrsz = get(0,'ScreenSize'); % Full screen
FnS_filter = 30; % How many WASP updates to plot
FnS = 82*FnS_filter;

map = loadmap;
set(map, 'OuterPosition', [1 1 scrsz(3)/2 scrsz(4)]);
hold on;
plot(r(1,1:FnS),r(2,1:FnS),'g-*')
%scatter(r(1,1:nS),r(2,1:nS),'g.') % Plot integrated
scatter(WASPloc(1,1:FnS_filter),WASPloc(2,1:FnS_filter),'rx') % Plot WASP observations
T = 1;
text(WASPloc(1,1:T:FnS_filter),WASPloc(2,1:T:FnS_filter),num2str(timestamps_WASP(1:T:FnS_filter)'))
legend('Anchor','Track','WASP observation');
hold off;

% Plot the WASP movements on the map (small version)