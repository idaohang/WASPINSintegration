%	WASPcheckData2.m
%
%	Mark Hedley
%	Created: 3 April 2009
%
% This function determines information about the frame structure and checks receive time consistency for errors.
% The following information from the data is returned in FrameInfo:
%	* SlotDuration			Duration (in seconds) of each slot
%	* SubSlotsPerSlot		Number of subslots in each slot
%	* SlotsPerSF			Number of slots in each superframe
%	* SubSlotDuration		Duration (in seconds) of each subslot
% The function uses the above info, along with the transmission slot to check that relative receive times are
% correct within a threshold given by SubSlotErrorMargin subslots. Errors are reported and bad data set as
% invalid.

function [FrameInfo SFdata] = WASPcheckData2(SFdata,NodeList,RefNodeIdx)

SubSlotErrorMargin = 100;

NumNodes = size(SFdata.TxSlot,2);
NumSF = SFdata.NumSF;

% Determine slot duration from data
SlotDurationList = [];
NumSF = SFdata.NumSF;
for SFidx=1:NumSF
	Vidx = find(squeeze(SFdata.RxValid(RefNodeIdx,:,SFidx)));
	if length(Vidx) < 3
		continue;
	end
	sRxTime = squeeze(SFdata.RxTime(:,:,SFidx));
	sSlot = squeeze(SFdata.TxSlot(SFidx,:));
	ListItem = (sRxTime(RefNodeIdx,Vidx(2:end)) - sRxTime(RefNodeIdx,Vidx(1)))./double((sSlot(Vidx(2:end)) - sSlot(Vidx(1))));
	SlotDurationList = [SlotDurationList ListItem];
	if length(SlotDurationList) > 1000
		break;
	end
end
SlotDuration = median(SlotDurationList);
SubSlotsPerSlot = round(SlotDuration/20.48e-6);
SlotDurationRound = SubSlotsPerSlot*20.48e-6;
disp(sprintf('INFO Slot length %.1f ms (%d sub-slots)',SlotDurationRound*1e3,SubSlotsPerSlot));

% determine SlotsPerSF
Vidx = find(SFdata.TxValid(:,RefNodeIdx));
SlotsPerSF = round(median(SFdata.TxTime(Vidx(2:end))-SFdata.TxTime(Vidx(1:(end-1))))/SlotDurationRound);
disp(sprintf('INFO %d slots per superframe',SlotsPerSF));

FrameInfo.SlotDuration = SlotDurationRound;
FrameInfo.SubSlotsPerSlot = SubSlotsPerSlot;
FrameInfo.SlotsPerSF = SlotsPerSF;
FrameInfo.SubSlotDuration = 20.48e-6;

% check that transmissions occur at the expected time.
% From the matrix of received times subtract the transmit time from the same node, thus each row is the
% relative receive time, with zero from self (i.e. each row is same but offset by relative transmit time). The
% expected value of this matrix is known from the transmit slot, so the error is the difference between these.
MaxTimeErr = SubSlotErrorMargin*FrameInfo.SubSlotDuration;
for SFidx=1:NumSF
	% extract SF data
	RxT = SFdata.RxTime(:,:,SFidx);
	RxV = SFdata.RxValid(:,:,SFidx);
	TxT = repmat(SFdata.TxTime(SFidx,:),NumNodes,1)';
	TxV = repmat(SFdata.TxValid(SFidx,:),NumNodes,1)';
	TxS = SFdata.TxSlot(SFidx,:);
	% determined matrix of expected time differences (probe delay between node pairs)
	ET = repmat(TxS,NumNodes,1);
	E2 = (ET - ET');
	% calculate error from expected time differences
	Tdiff = (RxT - TxT - double(E2)*SlotDurationRound).*RxV.*TxV;
	% report if any error too large - more than 2 subslots
	[v1 p1] = max(abs(Tdiff));
	[v2 p2] = max(v1);
	if v2 > MaxTimeErr
		% display error (only once per SF, not all bad values)
		disp(sprintf('WARNING: SFidx=%d Large error of %.2f ms (%.1f subslots) in receive time at node %d from %d', ...
			SFidx,Tdiff(p1(p2),p2)*1000, v2/FrameInfo.SubSlotDuration, NodeList(p1(p2)), NodeList(p2) ));
		% invalidate all data with large error
		GoodData = abs(Tdiff) <= MaxTimeErr;
		RxV = RxV.*GoodData;
		SFdata.RxValid(:,:,SFidx) = RxV;
	end
end
