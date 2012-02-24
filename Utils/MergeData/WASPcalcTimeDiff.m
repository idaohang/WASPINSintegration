function [TimeDiff TimeDiffValid] = WASPcalcTimeDiff(SFdata)
%	Written by Mark Hedley
%	Created 23/3/10
%
%	This function computes the relative time differences between the nodes (in seconds).

% compute difference between transmit/receive times for each beacon
NumNodes = size(SFdata.TxValid,2);
NumSF = SFdata.NumSF;
TD = zeros(NumNodes,NumNodes,NumSF);
TV = logical(zeros(NumNodes,NumNodes,NumSF));
for SFidx = 1:SFdata.NumSF
	
	RxT = squeeze(SFdata.RxTime(:,:,SFidx));
	RxV = logical(squeeze(SFdata.RxValid(:,:,SFidx)));
	TxT = squeeze(SFdata.TxTime(SFidx,:));
	TxV = logical(squeeze(SFdata.TxValid(SFidx,:)));
	
	for Nidx = 1:NumNodes
		if TxV(Nidx)==0
			continue;
		end
		
		Vrow = squeeze(RxV(:,Nidx));
		TD(Vrow,Nidx,SFidx) = RxT(Vrow,Nidx) - TxT(Nidx);
		TV(Vrow,Nidx,SFidx) = 1;
	end
end

% compute median over all valid values
MTD = zeros(NumNodes,NumNodes);
MTV = logical(zeros(NumNodes,NumNodes));
for N1=1:NumNodes
	for N2=1:NumNodes
		if sum(TV(N1,N2,:))<5
			continue;
		end
		MTD(N1,N2) = median(TD(N1,N2,TV(N1,N2,:)));
		MTV(N1,N2) = 1;
	end
end

% find column with most elements
TimeDiff = zeros(NumNodes,1);
TimeDiffValid = logical(zeros(NumNodes,1));
NodesToSearch = logical(ones(NumNodes,1));
OffsetVal = zeros(NumNodes,1);
Nelt = sum(MTV);
[v n] = max(Nelt);
if v>0
	TimeDiff = MTD(:,n);
	TimeDiffValid = MTV(:,n);
	NodesToSearch(n) = 0;
	OffsetVal(n) = 0;
else
	return;
end

% progressively find columns with offset for new node, and which has most overlap with those found
for NN=1:(NumNodes-1)
	% exit if offsets for all nodes found
	if min(TimeDiffValid)==1
		break;
	end

	% find columns with new values
	NewNodes = and(MTV,repmat(not(TimeDiffValid),1,NumNodes));
	NewCols = (sum(NewNodes) > 0);
	
	% find number of overlapping values
	ComNodes = and(MTV,repmat(TimeDiffValid,1,NumNodes));
	NumCommon = sum(ComNodes);
	
	% select column to match
	[v n] = max(NumCommon.*NewCols);
	if v==0
		return;
	else
		% determine offset from overlapping values
		AddedNodes = NewNodes(:,n);
		CommonNodes = ComNodes(:,n);
		OFS = median(MTD(CommonNodes,n) - TimeDiff(CommonNodes));
		TimeDiff(AddedNodes) = MTD(AddedNodes,n) - OFS;
		TimeDiffValid(AddedNodes) = 1;
		NodesToSearch(n) = 0;
		OffsetVal(n) = OFS;
	end
		
end

if min(TimeDiffValid) == 0
	disp(sprintf('WARNING: Time offset not found for all nodes'));
end

% determine an offset value for the remaining columns
NR = find(NodesToSearch)';
for n = NR
	CommonNodes = and(MTV(:,n),TimeDiffValid);
	if sum(CommonNodes) == 0
		continue;
	end
	OffsetVal(n) = median(MTD(CommonNodes,n) - TimeDiff(CommonNodes));
end

% create offset matrix
V1 = logical(MTV.*repmat(TimeDiffValid,1,NumNodes).*repmat(TimeDiffValid',NumNodes,1));
M1 = (MTD - repmat(OffsetVal',NumNodes,1)).*V1;

% for final time difference matrix take median of valid entries across rows of M1
TimeDiff2 = zeros(NumNodes,1);
VV = zeros(NumNodes,1);
for n=1:NumNodes
	v = V1(n,:);
	if length(v) > 0
		TimeDiff2(n) = median(M1(n,v));
		VV(n) = std(M1(n,v));
	end
end
TimeDiff = TimeDiff2 - min(TimeDiff2);
if max(abs(VV))>1e-4
	disp(sprintf('Inconsistency exceeds 0.1 ms'));
end
