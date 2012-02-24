% Perform maximum range gating
%
% Input:
%   tpos            - Location of the mobile node
%   anchorlocs      - Locations of the anchor nodes
%   meas            - Range measurement to each anchor.
%   max_range_error - Acceptable maximum range difference.
%
% Output:
%   anchorlocs      - Locations of the anchors satisfying the maximum range
%                     difference condition.
%   meas            - Rane measurements to corresponding anchors.
%
%
function [anchorlocs,meas,badAnchorIdx] = RangeGating(tpos, anchorlocs, meas, max_range_error)

badAnchorIdx = [];

num_sensors = length(meas);
prange = sqrt(sum((anchorlocs - tpos(:,ones(num_sensors,1))').^2,2));
nuk = meas - prange;

if any(abs(nuk)>max_range_error)
    badAnchorIdx = find(abs(nuk)>max_range_error);
    meas(badAnchorIdx) = [];
    anchorlocs(badAnchorIdx,:) = [];
end

% MAX_RANGE = 40;
% if any(abs(meas)>MAX_RANGE)
%     badAnchorIdx = find(abs(meas)>MAX_RANGE);
%     meas(badAnchorIdx) = [];
%     sensorPos(badAnchorIdx,:) = [];
% end

% badAnchorIdx = remnodes;
% meas(badAnchorIdx) = [];
% sensorPos(badAnchorIdx,:) = [];
