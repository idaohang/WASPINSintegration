% An implementation of the Grid Search algorithm
%
% Input:
%   r1          - Range measurements to the anchors
%   AnchorNodes - Locations of the anchor nodes
%   tpos        - Prior location of the mobile node (could be the predicted
%                 location)
%
% Output:
%   MobLocn_mean - Mean location calculated from the PDF on the grid
%   resid        - Residual of the location estiamte


function [MobLocn_mean, resid] = GridSearch(r1, AnchorNodes, tpos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these are only required for the get_locn_ml_mine.m function. it is
% defined here to speed up the script

if isempty(tpos)
    x_steps = 15:0.25:75;
    y_steps = -25:0.25:20;
else
    INT = 10;
    x_steps = tpos(1)-INT:0.1:tpos(1)+INT;
    y_steps = tpos(2)-INT:0.1:tpos(2)+INT;
end

total_probability = ones(length(y_steps),length(x_steps));
coordinate_gridx = repmat(x_steps,length(y_steps),1);
coordinate_gridy = repmat(y_steps',1,length(x_steps));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rlist = 1:length(r1);
forNonZeroProb1 = 0.01*ones(size(coordinate_gridx));

for r = rlist
    % for each valid range, construct the probability distribution on a grid
    this_range = r1(r);
    
    basexy = AnchorNodes(r,1:2);
    
    dist_grid = sqrt((coordinate_gridx - basexy(1)).^2 + (coordinate_gridy - basexy(2)).^2);
    error_grid = dist_grid - this_range;
    
    %Prob1 = ((error_grid < 0.2) & (error_grid > -2)) + forNonZeroProb1;
    Prob1 = ((error_grid < 0.5) & (error_grid > -4)) + forNonZeroProb1;
    
    total_probability = total_probability.*Prob1;
end

if all(all(total_probability == 0))
	MobLocn_mean = [];
    return;
end

% % maximum likelihood solution
% [y1,i1] = max(total_probability);
% [y2,xc1] = max(y1);
% yc1 = i1(xc1);
% ml_solution = [x_steps(xc1) y_steps(yc1)];
% MobLocn_ml = ml_solution;


% expected value solution
% Set total probability to zero where it is less than 1e-4 of the peak
probability_max = max(max(total_probability));
prob2 = (total_probability > (1e-4 * probability_max)) .* total_probability;
summed_probability = sum(sum(prob2));
Ex = sum(sum(coordinate_gridx .* prob2))/summed_probability;
Ey = sum(sum(coordinate_gridy .* prob2,2))/summed_probability;
MobLocn_mean = [Ex ; Ey];

resid = sum((r1 - sqrt(sum((AnchorNodes - repmat(MobLocn_mean',length(r1),1)).^2,2))).^2);

