function [ handle ] = loadmap( )
%loadmap Loads a map of the outdoor area from google maps
%   Also loads the positions of the base nodes from survey data


% PARAMETERS - The google map will be stretched to this.
min_x = -59.5;
max_x = 23;
min_y = -41;
% max_y is then determined by the scale of the map image

% Read map
map = imread('WASP_INS_Data/20111220/map.png','png');
max_y = min_y + size(map,1)/size(map,2)*(max_x-min_x);


% Read survey, flipping the x-axis
Survey = csvread('WASP_INS_Data/20111220/Survey_20Dec2011_with_fixed.csv');
Survey(9:12,:)=[]; % drop static surveys
anchors = num2str(Survey(:,1));
x = -Survey(:,2); %flip the x-axis
y = Survey(:,3);
z = Survey(:,4);


% Display, first flip the image upside down so it displays correctly
scrsz = get(0,'ScreenSize');
handle = figure('OuterPosition',scrsz);
imagesc([min_x max_x], [min_y max_y], flipdim(map,1));
axis equal;
axis([min_x max_x min_y max_y])
xlabel('x (m)');
ylabel('y (m)');
grid on;

hold on;
scatter(x,y,'filled')
text(x+1,y,anchors)

% now set the y-axis back to normal.
set(gca,'ydir','normal');
hold off



end

