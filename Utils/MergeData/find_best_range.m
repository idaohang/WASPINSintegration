% For static nodes finds the best range among the candidate ranges.
%
% Input:
%   RangeDataAll - Cell array of candidate ranges.
%   dist         - True Eucledian distance matrix between the nodes.
%
% Output:
%   BestRangeData - Ranges that match the true distances.
%
%

function BestRangeData = find_best_range(RangeDataAll,dist)

numNodes = size(dist,1);
numRanges = length(RangeDataAll);

for Nidx1 = 1:numNodes
    for Nidx2 = 1:numNodes
        BestRange = squeeze(RangeDataAll{1}(Nidx1,Nidx2,:));
        
        if all(BestRange==0)
            continue;
        end
        %BestRange(BestRange==0) = Inf;
        BestRangeErr = BestRange - dist(Nidx1,Nidx2);

        for Mtoa2 = 2:numRanges
            TempRange = squeeze(RangeDataAll{Mtoa2}(Nidx1,Nidx2,:));
            
            if all(TempRange==0)
                continue;
            end
            %TempRange(TempRange==0) = Inf;
            TempRangeErr = TempRange - dist(Nidx1,Nidx2);
            
            if any(abs(BestRangeErr)>abs(TempRangeErr))
                brlocn = find(abs(BestRangeErr)>abs(TempRangeErr));
                BestRange(brlocn) = TempRange(brlocn);
                BestRangeErr = BestRange - dist(Nidx1,Nidx2);
            end
        end

        RangeDataAll{1}(Nidx1,Nidx2,:) = BestRange;
    end
end

BestRangeData = RangeDataAll{1};
%plot_range_consistency(RangeDataAll{1}, RSSvalues, rxSNR, dist);
