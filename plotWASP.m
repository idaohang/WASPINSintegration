function [ handle ] = plotWASP( handle )

%plotWASP Plots a WASP csv data file on the map to just visualise the data
%   The handle is optional


% ask the user to specify a file
persistent pDir;
if isempty(pDir), pDir = '.'; end
[file, dirpath, fin] = uigetfile('*.csv', 'Pick a WASP data file', pDir);
if fin == 0 || isequal(file, 0) || isequal(dirpath, 0)
    return
end
pDir = dirpath;
filename = fullfile(pDir, file);

% Check if need to load the map
if nargin < 1
handle = loadmap();
end;

% now plot some data, again flipping the x-axis
figure(handle);
hold on;
% miss the first 4 observations which are frequently dodgy
wasp = csvread(filename,5,0);  
node220_x = -wasp(:,5);
node220_y = wasp(:,6);
node315_x = -wasp(:,12);
node315_y = wasp(:,13);
plot(node220_x,node220_y,'rs');
plot(node315_x,node315_y,'gs');
hold off;

% Calculate the distance between the nodes
error = sqrt(mean((node220_x - node315_x).^2 + (node220_y - node315_y).^2))


end

