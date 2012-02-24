% Converts MAC ID to Node IDs
% Input:
%   LogData - Data read using wtdread()
%
% Output:
%   LogData - Same input data with the MAC IDs replaced by Node IDs
%
% MacToID.csv is the file that contains the mapping.

function LogData = Convert_MAC_to_NodeID(LogData)

MacToIDFileName = 'MacToID.csv';

% sometimes all the elements of last row is zero, remove this
lengthLD = size(LogData,1);
remid = [];
for k1 = lengthLD:-1:1
    if all(LogData(k1,:) == 0)
        remid = [remid k1];
    else
        break;
    end
end
LogData(remid,:) = [];

% convert the mac ids to node ids
txnodes = unique(LogData(:,3));
rxnodes = unique(LogData(:,2));
mactoid = csvread(MacToIDFileName);
for k1 = 1:length(txnodes)
    nlocn = mactoid(:,1)==txnodes(k1);
    nid = mactoid(nlocn,2);
    rowids = LogData(:,3)==txnodes(k1);
    LogData(rowids,3) = nid;
end

for k1 = 1:length(rxnodes)
    nlocn = mactoid(:,1)==rxnodes(k1);
    nid = mactoid(nlocn,2);
    rowids = LogData(:,2)==rxnodes(k1);
    LogData(rowids,2) = nid;
end