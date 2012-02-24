% Remove data corresponding to certain nodes.
%
% Input:
%   NodesToRemove - IDs of the nodes, whose data to be removed.
%   Other input parameters are self-explanatory.
%
% Output:
%   Same as input, except the data corresponding to the specified nodes are
%   removed.
%
%

function [RangeDataAll, RSSvalues, rxSNR, TxTimes, TxTimesLocl, NodeLocn, NodeList, dist, ProbDataAll, MobileNodes] = ...
    remove_data(RangeDataAll, RSSvalues, rxSNR, TxTimes, TxTimesLocl, NodeLocn, NodeList, NodesToRemove, DIM, ProbDataAll, MobileNodes)

NodesToRemoveLocns = [];

for k1 = 1:length(NodesToRemove)
    
    nl = find(NodeList==NodesToRemove(k1));
    if ~isempty(nl)
        NodesToRemoveLocns = [NodesToRemoveLocns nl];
    end
end

RangeDataAll{1}(NodesToRemoveLocns,:,:) = [];
RangeDataAll{1}(:,NodesToRemoveLocns,:) = [];

RangeDataAll{2}(NodesToRemoveLocns,:,:) = [];
RangeDataAll{2}(:,NodesToRemoveLocns,:) = [];

RangeDataAll{3}(NodesToRemoveLocns,:,:) = [];
RangeDataAll{3}(:,NodesToRemoveLocns,:) = [];

RangeDataAll{4}(NodesToRemoveLocns,:,:) = [];
RangeDataAll{4}(:,NodesToRemoveLocns,:) = [];

RangeDataAll{5}(NodesToRemoveLocns,:,:) = [];
RangeDataAll{5}(:,NodesToRemoveLocns,:) = [];

RangeDataAll{6}(NodesToRemoveLocns,:,:) = [];
RangeDataAll{6}(:,NodesToRemoveLocns,:) = [];

RangeDataAll{7}(NodesToRemoveLocns,:,:) = [];
RangeDataAll{7}(:,NodesToRemoveLocns,:) = [];

RangeDataAll{8}(NodesToRemoveLocns,:,:) = [];
RangeDataAll{8}(:,NodesToRemoveLocns,:) = [];

RangeDataAll{9}(NodesToRemoveLocns,:,:) = [];
RangeDataAll{9}(:,NodesToRemoveLocns,:) = [];

ProbDataAll{1}(NodesToRemoveLocns,:,:) = [];
ProbDataAll{1}(:,NodesToRemoveLocns,:) = [];

ProbDataAll{2}(NodesToRemoveLocns,:,:) = [];
ProbDataAll{2}(:,NodesToRemoveLocns,:) = [];

ProbDataAll{3}(NodesToRemoveLocns,:,:) = [];
ProbDataAll{3}(:,NodesToRemoveLocns,:) = [];

ProbDataAll{4}(NodesToRemoveLocns,:,:) = [];
ProbDataAll{4}(:,NodesToRemoveLocns,:) = [];

ProbDataAll{5}(NodesToRemoveLocns,:,:) = [];
ProbDataAll{5}(:,NodesToRemoveLocns,:) = [];

ProbDataAll{6}(NodesToRemoveLocns,:,:) = [];
ProbDataAll{6}(:,NodesToRemoveLocns,:) = [];

ProbDataAll{7}(NodesToRemoveLocns,:,:) = [];
ProbDataAll{7}(:,NodesToRemoveLocns,:) = [];

ProbDataAll{8}(NodesToRemoveLocns,:,:) = [];
ProbDataAll{8}(:,NodesToRemoveLocns,:) = [];

ProbDataAll{9}(NodesToRemoveLocns,:,:) = [];
ProbDataAll{9}(:,NodesToRemoveLocns,:) = [];

RSSvalues(NodesToRemoveLocns,:,:) = [];
RSSvalues(:,NodesToRemoveLocns,:) = [];

rxSNR(NodesToRemoveLocns,:,:) = [];
rxSNR(:,NodesToRemoveLocns,:) = [];

TxTimes(NodesToRemoveLocns,:,:) = [];
TxTimesLocl(:,NodesToRemoveLocns) = [];

NodeLocn(NodesToRemoveLocns,:) = [];
NodeList(NodesToRemoveLocns,:) = [];

%MobileNodes(NodesToRemoveLocns) = [];

% calculate surveyed distances
NumNodes = size(NodeList,1);
if DIM == 2
    X = repmat(NodeLocn(:,2),1,NumNodes);
    Y = repmat(NodeLocn(:,3),1,NumNodes);
    
    dist = sqrt((X - X').^2 + (Y - Y').^2);
elseif DIM == 3
    X = repmat(NodeLocn(:,2),1,NumNodes);
    Y = repmat(NodeLocn(:,3),1,NumNodes);
    Z = repmat(NodeLocn(:,4),1,NumNodes);
    
    dist = sqrt((X - X').^2 + (Y - Y').^2 + (Z - Z').^2);
else
    error('Only 2 or 3 dimensional localization is practical.');
end

