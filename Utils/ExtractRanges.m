% Extract ranges that satisfy the minimum RSS and SNR conditions.
%
% Input:
%   LogFileNameBase - Filename of the data processed by MergeData.m
%                     function.
%   MobileNodes     - Node IDs of the mobile nodes.
%
% Output:
%   RangeData       - 3D array of range data (NumNodes*NumNodes*NumSF)   
%   NodeLocn        - Locations of all the nodes
%   MobileNodeLocn  - Location of the mobile node in the RangeData array.
%   
%

function [RangeData, NodeLocn, MobileNodeLocn] = ExtractRanges(LogFileName,MobileNodes)

MINRSS = -85;
MINSNR = 30;

% MINRSS = -90;
% MINSNR = 28;

load(LogFileName);

if ~isempty(MobileNodes)
    MobileNodeLocn = zeros(1,length(MobileNodes));
    for k1 = 1:length(MobileNodes)
        MobileNodeLocn(k1) = find(MobileNodes(k1)==NodeLocn(:,1));
    end
end
% if isempty(MobileNodes)
%     MobileNodes = NodeLocn(NodeLocn(:,end)==0,1);
% end
% 
% MobileNodeLocn = zeros(1,length(MobileNodes));
% for k1 = 1:length(MobileNodes)
%     MobileNodeLocn(k1) = find(MobileNodes(k1)==NodeLocn(:,1));
% end

% to use the best range data uncomment this line
% RangeData = BestRangeData;
numSF = size(RangeData,3);

% symmetrize the range data
for k1 = 1:numSF
    rd = RangeData(:,:,k1);
    
    RSSValid = (RSSvalues(:,:,k1) >= MINRSS) & (RSSvalues(:,:,k1) ~= 0);
    RSSValid = RSSValid .* RSSValid';
    rd = rd .* RSSValid;
    
    SNRValid = rxSNR(:,:,k1) >= MINSNR;
    SNRValid = SNRValid .* SNRValid';
    rd = rd .* SNRValid;
    
    rdvalid = rd >= 0.5;
    rd = rd .* rdvalid;
    
    RangeData(:,:,k1) = rd; %(rd + rd')/2;
end

if ~isempty(MobileNodes)
    % if there are any moile nodes that are not required to be tracked remove
    % them from RangeData
    MobileNodes1 = NodeLocn(NodeLocn(:,end)==0,1);
    MobileNodes1 = setdiff(MobileNodes1, MobileNodes);
    
    if ~isempty(MobileNodes1)
        MobileNodeLocn1 = zeros(1,length(MobileNodes1));
        for k1 = 1:length(MobileNodes1)
            MobileNodeLocn1(k1) = find(MobileNodes1(k1)==NodeLocn(:,1));
        end
        
        RangeData(MobileNodeLocn1,:,:) = [];
        RangeData(:,MobileNodeLocn1,:) = [];
        
        NodeLocn(MobileNodeLocn1,:) = [];
        
        MobileNodeLocn = zeros(1,length(MobileNodes));
        for k1 = 1:length(MobileNodes)
            MobileNodeLocn(k1) = find(MobileNodes(k1)==NodeLocn(:,1));
        end
    end
end

