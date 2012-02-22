function [ handle ] = loadfloorplan( )
%loadmap Loads a floorplan of the CSIRO Marsfield facility


% PARAMETERS - The floorplan will be stretched to this.
min_x = -35;
max_x = 39;
min_y = -35.3;
% max_y is then determined by the scale of the image

% Read map
map = imread('WASP_INS_Data/20120222/FloorPlan.png','png');
max_y = min_y + size(map,1)/size(map,2)*(max_x-min_x);

% Display, first flip the image upside down so it displays correctly
scrsz = get(0,'ScreenSize');
handle = figure('OuterPosition',scrsz);
imagesc([min_x max_x], [min_y max_y], flipdim(map,1));
axis equal;
axis([min_x max_x min_y max_y])
xlabel('x (m)');
ylabel('y (m)');

% now set the y-axis back to normal.
set(gca,'ydir','normal');
hold off



end

