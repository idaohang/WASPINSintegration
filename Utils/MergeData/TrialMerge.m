% Extract range measurements from raw data file. Can combine raw data from
% multiple collection nodes. It also extracts multiple candidate range
% measurements as well.
%
% Input:
%   Input parameters are self-explanatory.
%
% Output:
%   File with data is saved.
%
%

function TrialMerge(FileNameRxGainDelay,NodeDelayFileName,SurveyFileName,...
    FileNames,RangeDataFileName,NodesToRemove,DIM,REQUIRE_MAC_ID_TO_NODE_ID_CONVERSION,AIS_OR_ADHOC)

SpeedLight = 299792458*(1 - 0.00029);			% speed of light in vacuum, corrected for air at STP, in m/s
ExcessDelay = csvread(FileNameRxGainDelay)/2;	% additional receiver delay due to rx gain setting
NodeDelay = csvread(NodeDelayFileName);
LocnFixed = csvread(SurveyFileName);

if (AIS_OR_ADHOC==1) % AIS mode of data
    LocnFixed = LocnFixed';
    LocnFixed(LocnFixed(:,1)==0,:) = []; % remove all mobile nodes, they will be added from the data
    LocnFixed(:,1) = []; % remove first column
end

NumFiles = length(FileNames);

if NumFiles == 1
    LogData = wtdread(FileNames{1});
    if (REQUIRE_MAC_ID_TO_NODE_ID_CONVERSION)
        LogData = Convert_MAC_to_NodeID(LogData);
    end
    [comSFdata,uniqueNodes] = WASPgetSFdata2(LogData,ExcessDelay);
else
    % read files into memory
    for FI = 1:NumFiles
        % read and organise data from file
        % 	FileData{FI}.LogData = wtd(FileNames{FI});
        % 	[FileData{FI}.SFdata FileData{FI}.NodeList] = WASPgetSFdata2(FileData{FI}.LogData,ExcessDelay);
        LogData = wtdread(FileNames{FI});
        
        %     % sometimes all the elements of last row is zero, remove this
        %     if all(LogData(end,:) == 0)
        %         LogData(end,:) = [];
        %     end
        
        % replace the above lines with the following function
        if (REQUIRE_MAC_ID_TO_NODE_ID_CONVERSION)
            LogData = Convert_MAC_to_NodeID(LogData);
        end
        
        [FileData{FI}.SFdata FileData{FI}.NodeList] = WASPgetSFdata2(LogData,ExcessDelay);
        % find equivalent transmit time (in slot zero) for each transmitter
        NodeList = FileData{FI}.NodeList;
        NumNodes = length(NodeList);
        NumSF = FileData{FI}.SFdata.NumSF;
        TxTime = FileData{FI}.SFdata.TxTime;
        TxValid = FileData{FI}.SFdata.TxValid;
        TxSlot = FileData{FI}.SFdata.TxSlot;
        SeqNum = FileData{FI}.SFdata.SeqNum;
        SyncNodeInfo = [];
        SFduration = [];
        for NI = 1:NumNodes
            % find first valid transmit time, if any
            Vidx = find(TxTime(:,NI));
            if isempty(Vidx)
                continue;		% no valid tx time for node, move onto next node
            end
            % record node data
            SyncNodeInfo = [SyncNodeInfo ; [double(NodeList(NI)) double(TxTime(Vidx(1),NI))...
                double(TxSlot(Vidx(1),NI)) double(Vidx(1)) double(SeqNum(Vidx(1))) ] ];
            % update data on SF duration if consecutive tx times
            if (length(Vidx) >=2) && (SeqNum(Vidx(2)) == SeqNum(Vidx(1))+1) && (TxSlot(Vidx(2),NI) == TxSlot(Vidx(1),NI))
                SFduration = [SFduration ; (TxTime(Vidx(2),NI)-TxTime(Vidx(1),NI))];
            end
        end
        
        % estimate duration superframe
        FileData{FI}.PeriodSF = median(SFduration);
        
        % estimate TxTime for each node in first SF
        FirstSlotTxTime = SyncNodeInfo(:,2) - (SyncNodeInfo(:,4)-1) * FileData{FI}.PeriodSF;
        SyncNodeInfo = [SyncNodeInfo FirstSlotTxTime];
        
        FileData{FI}.SyncNodeInfo = SyncNodeInfo;
        FileData{FI}.SFduration = SFduration;
    end
    
    % match files
    FileSFoffset = zeros(NumFiles,NumFiles);
    FileSFoffsetValid = zeros(NumFiles,NumFiles);
    A = [];
    b = [];
    
    % try all combinations of files
    for F1 = 1:NumFiles
        for F2 = (F1+1):NumFiles
            % - find common nodes
            CommonNodes = intersect(FileData{F1}.SyncNodeInfo(:,1),FileData{F2}.SyncNodeInfo(:,1));
            if isempty(CommonNodes)
                continue;
            end
            % - find index for first common node into each file
            V1 = find(FileData{F1}.SyncNodeInfo(:,1) == CommonNodes(1));
            V1 = V1(1);
            V2 = find(FileData{F2}.SyncNodeInfo(:,1) == CommonNodes(1));
            V2 = V2(1);
            FileSFoffset(F1,F2) = (FileData{F2}.SyncNodeInfo(V2,6) - FileData{F1}.SyncNodeInfo(V1,6))/FileData{FI}.PeriodSF;
            FileSFoffsetValid(F1,F2) = 1;
            
            % form consistency matrix
            RD = zeros(1,NumFiles);
            RD(F1) = 1;
            RD(F2) = -1;
            A = [A ; RD];
            b = [b ; round(FileSFoffset(F1,F2))];
        end
    end
    
    FileSFoffset = round(FileSFoffset);
    
    % construct and solve matrix for relative time
    A2 = A(:,2:end);
    x = [ 0 ; inv(A2'*A2)*A2'*b];
    x2 = x + 1 - min(x);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    uniqueNodes = [];
    noSFInFile = zeros(1,NumFiles);
    for FI = 1:NumFiles
        noSFInFile(FI) = FileData{FI}.SFdata.NumSF;
        uniqueNodes = [uniqueNodes ; FileData{FI}.NodeList];
    end
    %maxNoSF = max(noSFInFile);
    maxNoSF = min(noSFInFile); % take the minimum number of superframes rather than the maximum number
    uniqueNodes = unique(uniqueNodes);
    numUniqueNodes = length(uniqueNodes);
    
    nodeMap = zeros(numUniqueNodes,NumFiles);
    for FI = 1:NumFiles
        fNodeList = FileData{FI}.NodeList;
        fNumNodes = length(fNodeList);
        for fiNidx = 1:fNumNodes
            cNidx = fNodeList(fiNidx);
            cNidxLocn = find(uniqueNodes==cNidx);
            if ~isempty(cNidxLocn)
                nodeMap(cNidxLocn,FI) = 1;
            end
        end
    end
    
    comSFdata = [];
    
    tNumSF = 0;
    curSeqNum = zeros(1,NumFiles);
    
    % allocate space
    comSFdata.NumSF = 0;
    comSFdata.SeqNum = zeros(maxNoSF,1);
    
    comSFdata.TxValid = zeros(maxNoSF,numUniqueNodes);
    comSFdata.TxTime = zeros(maxNoSF,numUniqueNodes);
    comSFdata.TxSlot = zeros(maxNoSF,numUniqueNodes);
    comSFdata.TxGain = zeros(maxNoSF,numUniqueNodes);
    comSFdata.Vbat = zeros(maxNoSF,numUniqueNodes);
    
    comSFdata.RxValid = zeros(numUniqueNodes,numUniqueNodes,maxNoSF);
    comSFdata.RxTime = zeros(numUniqueNodes,numUniqueNodes,maxNoSF);
    comSFdata.RxTime2 = zeros(numUniqueNodes,numUniqueNodes,maxNoSF);
    comSFdata.RxTimeMax = zeros(numUniqueNodes,numUniqueNodes,maxNoSF);
    comSFdata.RxAmp = zeros(numUniqueNodes,numUniqueNodes,maxNoSF);
    comSFdata.RxAmp2 = zeros(numUniqueNodes,numUniqueNodes,maxNoSF);
    comSFdata.RxAmpMax = zeros(numUniqueNodes,numUniqueNodes,maxNoSF);
    comSFdata.RxGain = zeros(numUniqueNodes,numUniqueNodes,maxNoSF);
    comSFdata.RxSNR = zeros(numUniqueNodes,numUniqueNodes,maxNoSF);
    
    skipFile = zeros(1,NumFiles);
    
    for SFidx = 1:maxNoSF
        
        cSFidx = SFidx - 1 + x2; %
        % check whether this SF exists in all files
        if (any(cSFidx > maxNoSF))
            break;
        end
        
        tNumSF = tNumSF + 1;
        for FI = 1:NumFiles
            curSeqNum(FI) = FileData{FI}.SFdata.SeqNum(cSFidx(FI));
        end
        
        if any(curSeqNum - curSeqNum(1))
            fprintf('Current sequence numbers in different files are not the same.\n');
            locn = find(curSeqNum - min(curSeqNum) ~= 0);
            x2(locn) = x2(locn) - 1;
            skipFile(locn) = 1;
            %continue;
        end
        
        comSFdata.NumSF = tNumSF;
        comSFdata.SeqNum(tNumSF,1) = curSeqNum(1);
        
        % transmit side of things
        txv = zeros(length(uniqueNodes),NumFiles);
        txt = zeros(length(uniqueNodes),NumFiles);
        txs = zeros(length(uniqueNodes),NumFiles);
        txg = zeros(length(uniqueNodes),NumFiles);
        txvbat = zeros(length(uniqueNodes),NumFiles);
        
        for FI = 1:NumFiles
            if ~skipFile(FI)
                txv(nodeMap(:,FI)==1,FI) = FileData{FI}.SFdata.TxValid(cSFidx(FI),:);
                txt(nodeMap(:,FI)==1,FI) = FileData{FI}.SFdata.TxTime(cSFidx(FI),:);
                txs(nodeMap(:,FI)==1,FI) = FileData{FI}.SFdata.TxSlot(cSFidx(FI),:);
                txg(nodeMap(:,FI)==1,FI) = FileData{FI}.SFdata.TxGain(cSFidx(FI),:);
                txvbat(nodeMap(:,FI)==1,FI) = FileData{FI}.SFdata.Vbat(cSFidx(FI),:);
            end
        end
        
        comSFdata.TxValid(tNumSF,:) = any(txv')';
        
        txt1 = (txt - repmat(max(txt,[],2),1,NumFiles)) .* txv;
        if (any(any(txt1)))
            fprintf('Not all recorded transmit times of at least one of the node in the current superframe are not equal.\n');
        else
            comSFdata.TxTime(tNumSF,:) = max(txt,[],2);
        end
        txs1 = (txs - repmat(max(txs,[],2),1,NumFiles)) .* txv;
        if (any(any(txs1)))
            fprintf('Not all recorded tx slot of at least one of the node in the current superframe are not equal.\n');
        else
            comSFdata.TxSlot(tNumSF,:) = max(txs,[],2);
        end
        txg1 = (txg - repmat(max(txg,[],2),1,NumFiles)) .* txv;
        if (any(any(txg1)))
            %fprintf('Not all recorded transmit gains of at least one of the node in the current superframe are not equal.\n');
        else
            comSFdata.TxGain(tNumSF,:) = max(txg,[],2);
        end
        txvbat1 = (txvbat - repmat(max(txvbat,[],2),1,NumFiles)) .* txv;
        if (any(any(txvbat1)))
            %fprintf('Not all recorded battery voltages of at least one of the node in the current superframe are not equal.\n');
        else
            comSFdata.Vbat(tNumSF,:) = max(txvbat,[],2);
        end
        
        % receive side of things
        rxv = zeros(length(uniqueNodes),length(uniqueNodes));
        rxt = zeros(length(uniqueNodes),length(uniqueNodes));
        rxa = zeros(length(uniqueNodes),length(uniqueNodes));
        rxt_2 = zeros(length(uniqueNodes),length(uniqueNodes));
        rxa_2 = zeros(length(uniqueNodes),length(uniqueNodes));
        rxt_m = zeros(length(uniqueNodes),length(uniqueNodes));
        rxa_m = zeros(length(uniqueNodes),length(uniqueNodes));
        rxg = zeros(length(uniqueNodes),length(uniqueNodes));
        rxsnr = zeros(length(uniqueNodes),length(uniqueNodes));
        
        startFid = 1;
        for FI = 1:NumFiles
            if skipFile(FI)
                startFid = FI + 1;
                continue;
            end
            rxv(nodeMap(:,FI)==1,nodeMap(:,FI)==1) = FileData{FI}.SFdata.RxValid(:,:,cSFidx(FI));
            rxt(nodeMap(:,FI)==1,nodeMap(:,FI)==1) = FileData{FI}.SFdata.RxTime(:,:,cSFidx(FI));
            rxa(nodeMap(:,FI)==1,nodeMap(:,FI)==1) = FileData{FI}.SFdata.RxAmp(:,:,cSFidx(FI));
            rxt_2(nodeMap(:,FI)==1,nodeMap(:,FI)==1) = FileData{FI}.SFdata.RxTime2(:,:,cSFidx(FI));
            rxa_2(nodeMap(:,FI)==1,nodeMap(:,FI)==1) = FileData{FI}.SFdata.RxAmp2(:,:,cSFidx(FI));
            rxt_m(nodeMap(:,FI)==1,nodeMap(:,FI)==1) = FileData{FI}.SFdata.RxTimeMax(:,:,cSFidx(FI));
            rxa_m(nodeMap(:,FI)==1,nodeMap(:,FI)==1) = FileData{FI}.SFdata.RxAmpMax(:,:,cSFidx(FI));
            rxg(nodeMap(:,FI)==1,nodeMap(:,FI)==1) = FileData{FI}.SFdata.RxGain(:,:,cSFidx(FI));
            rxsnr(nodeMap(:,FI)==1,nodeMap(:,FI)==1) = FileData{FI}.SFdata.RxSNR(:,:,cSFidx(FI));
            
            if (FI == startFid)
                rxv0 = rxv;
                rxt0 = rxt;
                rxa0 = rxa;
                rxt_20 = rxt_2;
                rxa_20 = rxa_2;
                rxt_m0 = rxt_m;
                rxa_m0 = rxa_m;
                rxg0 = rxg;
                rxsnr0 = rxsnr;
            else
                rxt2 = rxt0 - rxt;
                if any(any(rxt2 .* rxv .* rxv0))
                    fprintf('Not all recorded received times of at least one of the node in the current superframe are not equal.\n');
                end
                rxt2 = rxt0 + rxt - (rxt .* rxv .*rxv0);
                
                rxa2 = rxa0 - rxa;
                if any(any(rxa2 .* rxv .* rxv0))
                    fprintf('Not all recorded received amplitude of at least one of the node in the current superframe are not equal.\n');
                end
                rxa2 = rxa0 + rxa - (rxa .* rxv .*rxv0);
                
                rxt_22 = rxt_20 - rxt_2;
                if any(any(rxt_22 .* rxv .* rxv0))
                    fprintf('Not all recorded received times of at least one of the node in the current superframe are not equal.\n');
                end
                rxt_22 = rxt_20 + rxt_2 - (rxt_2 .* rxv .*rxv0);
                
                rxa_22 = rxa_20 - rxa_2;
                if any(any(rxa_22 .* rxv .* rxv0))
                    fprintf('Not all recorded received times of at least one of the node in the current superframe are not equal.\n');
                end
                rxa_22 = rxa_20 + rxa_2 - (rxa_2 .* rxv .*rxv0);
                
                rxt_m2 = rxt_m0 - rxt_m;
                if any(any(rxt_m2 .* rxv .* rxv0))
                    fprintf('Not all recorded received times of at least one of the node in the current superframe are not equal.\n');
                end
                rxt_m2 = rxt_m0 + rxt_m - (rxt_m .* rxv .*rxv0);
                
                rxa_m2 = rxa_m0 - rxa_m;
                if any(any(rxa_m2 .* rxv .* rxv0))
                    fprintf('Not all recorded received times of at least one of the node in the current superframe are not equal.\n');
                end
                rxa_m2 = rxa_m0 + rxa_m - (rxa_m .* rxv .*rxv0);
                
                rxg2 = rxg0 - rxg;
                if any(any(rxg2 .* rxv .* rxv0))
                    fprintf('Not all recorded received times of at least one of the node in the current superframe are not equal.\n');
                end
                rxg2 = rxg0 + rxg - (rxg .* rxv .*rxv0);
                
                rxsnr2 = rxsnr0 - rxsnr;
                if any(any(rxsnr2 .* rxv .* rxv0))
                    fprintf('Not all recorded received times of at least one of the node in the current superframe are not equal.\n');
                end
                rxsnr2 = rxsnr0 + rxsnr - (rxsnr .* rxv .*rxv0);
                
                rxv0 = or(rxv,rxv0);
                rxt0 = rxt2;
                rxa0 = rxa2;
                rxt_20 = rxt_22;
                rxa_20 = rxa_22;
                rxt_m0 = rxt_m2;
                rxa_m0 = rxa_m2;
                rxg0 = rxg2;
                rxsnr0 = rxsnr2;
            end
        end
        
        skipFile = zeros(1,NumFiles);
        
        comSFdata.RxValid(:,:,tNumSF) = rxv0;
        comSFdata.RxTime(:,:,tNumSF) = rxt0;
        comSFdata.RxTime2(:,:,tNumSF) = rxt_20;
        comSFdata.RxTimeMax(:,:,tNumSF) = rxt_m0;
        comSFdata.RxAmp(:,:,tNumSF) = rxa0;
        comSFdata.RxAmp2(:,:,tNumSF) = rxa_20;
        comSFdata.RxAmpMax(:,:,tNumSF) = rxa_m0;
        comSFdata.RxGain(:,:,tNumSF) = rxg0;
        comSFdata.RxSNR(:,:,tNumSF) = rxsnr0;
    end
    
    comSFdata.SeqNum(tNumSF+1:end) = [];
    comSFdata.TxValid(tNumSF+1:end,:) = [];
    comSFdata.TxTime(tNumSF+1:end,:) = [];
    comSFdata.TxSlot(tNumSF+1:end,:) = [];
    comSFdata.TxGain(tNumSF+1:end,:) = [];
    comSFdata.Vbat(tNumSF+1:end,:) = [];
    
    comSFdata.RxValid(:,:,tNumSF+1:end) = [];
    comSFdata.RxTime(:,:,tNumSF+1:end) = [];
    comSFdata.RxTime2(:,:,tNumSF+1:end) = [];
    comSFdata.RxTimeMax(:,:,tNumSF+1:end) = [];
    comSFdata.RxAmp(:,:,tNumSF+1:end) = [];
    comSFdata.RxAmp2(:,:,tNumSF+1:end) = [];
    comSFdata.RxAmpMax(:,:,tNumSF+1:end) = [];
    comSFdata.RxGain(:,:,tNumSF+1:end) = [];
    comSFdata.RxSNR(:,:,tNumSF+1:end) = [];  
end


NodeList = uniqueNodes;
NumNodes = length(NodeList);
% arbitrarily pick node as reference node
RefNodeIdx = 1;

% check data
%[FrameInfo SFdata] = WASPcheckData2(comSFdata,NodeList,RefNodeIdx);
[FrameInfo,SFdata] = WASPcheckData2(comSFdata,NodeList,RefNodeIdx);

% calculate time offsets of each node and obtain the synchronized
% localization txtime of each node
[TimeDiff TimeDiffValid] = WASPcalcTimeDiff(SFdata);
if (any(TimeDiffValid==0))
    fprintf('WARNING: TxOffset calculated is not valid for at least one of the nodes.\n');
end
TxTimes = SFdata.TxTime' - repmat(TimeDiff,1,SFdata.NumSF);

% calculate round trip delay for node pairs
[RTDdata RSSvalues] = WASPcalcRTD2(SFdata,1);

% % update delay, check same nodes for all files
ncols = size(LocnFixed,2);
NodeLocn = zeros(NumNodes,ncols+1);
MobileNodes = [];
for NI = 1:NumNodes
    nodeidx = find(LocnFixed(:,1)==NodeList(NI),1);
    
    if ~isempty(nodeidx)
        NodeLocn(NI,1:ncols) = LocnFixed(LocnFixed(:,1)==NodeList(NI),:);
        NodeLocn(NI,ncols+1) = NodeDelay(NodeDelay(:,1)==NodeList(NI),2);
    else
        fprintf('A mobile node (ID %d) is found, setting its position to zero in the node list.\n', NodeList(NI));
        NodeLocn(NI,1:ncols) = [NodeList(NI) zeros(1,ncols-1)];
        NodeLocn(NI,ncols+1) = NodeDelay(NodeDelay(:,1)==NodeList(NI),2);
        
        MobileNodes = [MobileNodes NodeList(NI)];
    end
    %     VI = find(LocnFixed(:,1)==NodeList(NI));
    %     NodeLocn(NI,1:3) = LocnFixed(VI,:);
    %     VI = find(NodeDelay(:,1)==NodeList(NI));
    %     NodeLocn(NI,5) = NodeDelay(VI,2);
end

% compute range
RangeData  = RTDdata.RTD;		% range using first TOA from both nodes
RangeData2 = squeeze(RTDdata.RTD5(:,:,1,:));
RangeData3 = squeeze(RTDdata.RTD5(:,:,2,:));
RangeData4 = squeeze(RTDdata.RTD5(:,:,3,:));
RangeData5 = squeeze(RTDdata.RTD5(:,:,4,:));
RangeData6 = squeeze(RTDdata.RTD5(:,:,5,:));
RangeData7 = squeeze(RTDdata.RTD5(:,:,6,:));
RangeData8 = squeeze(RTDdata.RTD5(:,:,7,:));
RangeData9 = squeeze(RTDdata.RTD5(:,:,8,:));

ProbData  = RTDdata.Prob;		% probability each range is a true range
ProbData2 = squeeze(RTDdata.Prob5(:,:,1,:));
ProbData3 = squeeze(RTDdata.Prob5(:,:,2,:));
ProbData4 = squeeze(RTDdata.Prob5(:,:,3,:));
ProbData5 = squeeze(RTDdata.Prob5(:,:,4,:));
ProbData6 = squeeze(RTDdata.Prob5(:,:,5,:));
ProbData7 = squeeze(RTDdata.Prob5(:,:,6,:));
ProbData8 = squeeze(RTDdata.Prob5(:,:,7,:));
ProbData9 = squeeze(RTDdata.Prob5(:,:,8,:));

for Nidx = 1:NumNodes
    RangeData(:,Nidx,:) = RangeData(:,Nidx,:) - NodeLocn(Nidx,end);
    RangeData(Nidx,:,:) = RangeData(Nidx,:,:) - NodeLocn(Nidx,end);    
    RangeData2(:,Nidx,:) = RangeData2(:,Nidx,:) - NodeLocn(Nidx,end);
    RangeData2(Nidx,:,:) = RangeData2(Nidx,:,:) - NodeLocn(Nidx,end);
    RangeData3(:,Nidx,:) = RangeData3(:,Nidx,:) - NodeLocn(Nidx,end);
    RangeData3(Nidx,:,:) = RangeData3(Nidx,:,:) - NodeLocn(Nidx,end);
    RangeData4(:,Nidx,:) = RangeData4(:,Nidx,:) - NodeLocn(Nidx,end);
    RangeData4(Nidx,:,:) = RangeData4(Nidx,:,:) - NodeLocn(Nidx,end);
    RangeData5(:,Nidx,:) = RangeData5(:,Nidx,:) - NodeLocn(Nidx,end);
    RangeData5(Nidx,:,:) = RangeData5(Nidx,:,:) - NodeLocn(Nidx,end);
    RangeData6(:,Nidx,:) = RangeData6(:,Nidx,:) - NodeLocn(Nidx,end);
    RangeData6(Nidx,:,:) = RangeData6(Nidx,:,:) - NodeLocn(Nidx,end);
    RangeData7(:,Nidx,:) = RangeData7(:,Nidx,:) - NodeLocn(Nidx,end);
    RangeData7(Nidx,:,:) = RangeData7(Nidx,:,:) - NodeLocn(Nidx,end);
    RangeData8(:,Nidx,:) = RangeData8(:,Nidx,:) - NodeLocn(Nidx,end);
    RangeData8(Nidx,:,:) = RangeData8(Nidx,:,:) - NodeLocn(Nidx,end);
    RangeData9(:,Nidx,:) = RangeData9(:,Nidx,:) - NodeLocn(Nidx,end);
    RangeData9(Nidx,:,:) = RangeData9(Nidx,:,:) - NodeLocn(Nidx,end);
end
RangeData  = (RangeData  * SpeedLight/2).*RTDdata.RTDvalid;
RangeData2 = (RangeData2 * SpeedLight/2).*RTDdata.RTDvalid;
RangeData3 = (RangeData3 * SpeedLight/2).*RTDdata.RTDvalid;
RangeData4 = (RangeData4 * SpeedLight/2).*RTDdata.RTDvalid;
RangeData5 = (RangeData5 * SpeedLight/2).*RTDdata.RTDvalid;
RangeData6 = (RangeData6 * SpeedLight/2).*RTDdata.RTDvalid;
RangeData7 = (RangeData7 * SpeedLight/2).*RTDdata.RTDvalid;
RangeData8 = (RangeData8 * SpeedLight/2).*RTDdata.RTDvalid;
RangeData9 = (RangeData9 * SpeedLight/2).*RTDdata.RTDvalid;

rxSNR = double(SFdata.RxSNR) .* RTDdata.RTDvalid;
RSSvalues = RSSvalues .* RTDdata.RTDvalid;

numSF = size(RangeData,3);
remid = zeros(1,numSF);

for SFidx = 1:numSF
    rd = RangeData(:,:,SFidx);
    
    if all(rd(:)==0)
        remid(1,SFidx) = SFidx;
    end
end
remid(remid==0) = [];
RangeData(:,:,remid) = [];
RangeData2(:,:,remid) = [];
RangeData3(:,:,remid) = [];
RangeData4(:,:,remid) = [];
RangeData5(:,:,remid) = [];
RangeData6(:,:,remid) = [];
RangeData7(:,:,remid) = [];
RangeData8(:,:,remid) = [];
RangeData9(:,:,remid) = [];

ProbData(:,:,remid) = [];
ProbData2(:,:,remid) = [];
ProbData3(:,:,remid) = [];
ProbData4(:,:,remid) = [];
ProbData5(:,:,remid) = [];
ProbData6(:,:,remid) = [];
ProbData7(:,:,remid) = [];
ProbData8(:,:,remid) = [];
ProbData9(:,:,remid) = [];

RSSvalues(:,:,remid) = [];
rxSNR(:,:,remid) = [];
TxTimes(:,remid) = [];

% synchronized times and local clock transmission times
TxTimesSync = TxTimes;
TxTimesLocl = SFdata.TxTime;

TxTimesLocl(remid,:) = [];

RangeDataAll{1} = RangeData;
RangeDataAll{2} = RangeData2;
RangeDataAll{3} = RangeData3;
RangeDataAll{4} = RangeData4;
RangeDataAll{5} = RangeData5;
RangeDataAll{6} = RangeData6;
RangeDataAll{7} = RangeData7;
RangeDataAll{8} = RangeData8;
RangeDataAll{9} = RangeData9;

ProbDataAll{1} = ProbData;
ProbDataAll{2} = ProbData2;
ProbDataAll{3} = ProbData3;
ProbDataAll{4} = ProbData4;
ProbDataAll{5} = ProbData5;
ProbDataAll{6} = ProbData6;
ProbDataAll{7} = ProbData7;
ProbDataAll{8} = ProbData8;
ProbDataAll{9} = ProbData9;

[RangeDataAll, RSSvalues, rxSNR, TxTimesSync, TxTimesLocl, NodeLocn, NodeList, dist, ProbDataAll, MobileNodes] = ...
    remove_data(RangeDataAll, RSSvalues, rxSNR, TxTimesSync, TxTimesLocl, NodeLocn, NodeList, NodesToRemove, DIM, ProbDataAll, MobileNodes);
RangeData = RangeDataAll{1};

MobileNodeLocn = [];
for k1 = 1:length(MobileNodes)
    nl = find(MobileNodes(k1)==NodeList);
    if ~isempty(nl)
        MobileNodeLocn = [MobileNodeLocn nl];
    end
end

BestRangeData = []; %find_best_range(RangeDataAll,dist); % only valid for static nodes

%plot_range_consistency(RangeData, RSSvalues, rxSNR, dist);

% last column is node time delay not required anymore, remove it
NodeLocn(:,end) = [];

save(RangeDataFileName, 'RangeData','RangeDataAll', 'RSSvalues', 'rxSNR', ...
    'TxTimesSync', 'TxTimesLocl', 'NodeLocn', 'NodeList', 'dist', 'BestRangeData', 'ProbDataAll','MobileNodeLocn');


