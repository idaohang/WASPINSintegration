%	WASPcalcRTD2.m
%
%	Mark Hedley
%	Created: 3 April 2009
%	Modified: 19 May 2009
%	Modified: 19 Jan 2010 to support changes in data type in SFdata, and calculate all RTD values
%
% This function takes the data returned by WASPgetSFdata and computes the round trip delay between node pairs.
% The frequency difference between node pairs is computed as an intermediate variable.
% Output is:
%	* RTDdata.FreqDiff			frequency difference (fractional) between node pairs
%	* RTDdata.FreqDiffValid		binary flag, 1 for valid values in RTDdata.FreqDiff
%	* RTDdata.RTD				round trip delay between node pair (in seconds), corrected for probe delay,
%								using primary TOA values
%	* RTDdata.RTD5				RTD values using other TOA values, only returned when CalcAllRTD is non-zero
%								Order: 1:1-2 2:1-Max 3:2-2 4:2-Max 5:Max-Max
%	* RTDdata.RTDvalid			binary flag, 1 for valid values in RTDdata.RTD
%	* RSSvalues					received signal strength (in dBm), only valid for SFdata.RxValid

function [RTDdata RSSvalues] = WASPcalcRTD2(SFdata,CalcAllRTD)


prob_table1 = 2*[30 42 60 85 120 169 200 220]; % SNRs (ratio*16) 
prob_table2 = [96 160 224 288 352 416 500 580]; % times (in units of 5/32ns)
good_prob = ...
[0.17 0.15 0.12 0.10 0.08 0.01 0.01 0.01;...
0.17 0.15 0.18 0.05 0.01 0.01 0.01 0.01;...
0.26 0.38 0.47 0.48 0.1 0.01 0.01 0.01;...
0.70 0.70 0.95 0.90 0.1 0.01 0.01 0.01;...
0.95 0.93 0.93 0.95 0.2 0.01 0.01 0.01;...
0.99 0.99 0.99 0.95 0.5 0.1 0.1 0.01;...
0.99 0.99 0.99 0.98 0.95 0.4 0.3 0.3;...
0.99 0.99 0.99 0.99 0.98 0.49 0.4 0.4];


FreqDiffErrThres = 1e-6;			% threshold for error in frequency difference

NumNodes = size(SFdata.TxSlot,2);
NumSF = SFdata.NumSF;

% initialise frequency difference and RTD data
RTDdata.FreqDiff = zeros(NumNodes,NumNodes,NumSF);
RTDdata.FreqDiffValid = false(NumNodes,NumNodes,NumSF);
RTDdata.RTD = zeros(NumNodes,NumNodes,NumSF);
RTDdata.Prob = zeros(NumNodes,NumNodes,NumSF);

RTDdata.RTDvalid = false(NumNodes,NumNodes,NumSF);
if CalcAllRTD~=0
	%RTDdata.RTD5 = zeros(NumNodes,NumNodes,5,NumSF);
    RTDdata.RTD5 = zeros(NumNodes,NumNodes,8,NumSF);
    RTDdata.Prob5 = zeros(NumNodes,NumNodes,8,NumSF);
end

% calculate received signal strength in dBm (roughly)
RSSvalues = -((SFdata.RxGain >= 96)*34.5 + 2*double(bitand(SFdata.RxGain,31)) + 18).*SFdata.RxValid;

% loop through SF computing RTD, invalid for first SF so skipped
for SFidx=2:NumSF
	
	% extract rx time and flag for current and previous SF
	CurrentRxValid = squeeze(SFdata.RxValid(:,:,SFidx));
	PrevRxValid = squeeze(SFdata.RxValid(:,:,SFidx-1));
	CurrentRxTime = squeeze(SFdata.RxTime(:,:,SFidx));
	PrevRxTime = squeeze(SFdata.RxTime(:,:,SFidx-1));

	% determine probe delay
	%ET = repmat(SFdata.TxSlot(SFidx,:),NumNodes,1);		% don't use slot for variable payload length
	%ProbeDelay = (ET - ET')*FrameInfo.SlotDuration;
	PDAi = repmat(SFdata.TxTime(SFidx,:),NumNodes,1)+squeeze(SFdata.RxTime(:,:,SFidx));
	ProbeDelay = (PDAi-PDAi')/2;	% Col Node Tx Time - Row Node Tx Time

	% determine frequency difference
	FDValid = and(CurrentRxValid,PrevRxValid);
	TxDiff = SFdata.TxTime(SFidx,:) - SFdata.TxTime(SFidx-1,:);
	TxDiff = TxDiff + (TxDiff == 0);	% prevent divide by zero
	TxDiff2 = repmat(TxDiff,NumNodes,1);
	TxTimeArray = repmat(SFdata.TxTime(SFidx,:),NumNodes,1);
	CurrentTT = CurrentRxTime - TxTimeArray;
	PrevTT = PrevRxTime - repmat(SFdata.TxTime(SFidx-1,:),NumNodes,1);
	CurrentFreqDiff = (CurrentTT-PrevTT)./TxDiff2.*FDValid;
	
	% display unexpectedly large frequency differences, but only if signal strength >= -85 dBm
	if max(max(abs(CurrentFreqDiff))) > FreqDiffErrThres
		RSS = RSSvalues(:,:,SFidx);
		MFD = CurrentFreqDiff.*(RSS > -85);
		[a b] = max(abs(MFD));
		[c d] = max(a);
		if c > FreqDiffErrThres
			disp(sprintf('Big FreqDiff error at %d,%d, RSS=%.0f,%.0f',d,b(d),RSS(d,b(d)),RSS(b(d),d)));
		end
	end
	
	FDValid = and(FDValid,(abs(CurrentFreqDiff) < FreqDiffErrThres));	% reject large values
	RTDdata.FreqDiff(:,:,SFidx) = CurrentFreqDiff;
	RTDdata.FreqDiffValid(:,:,SFidx) = FDValid;
	
	% PHIL to take a close look at what is happening you can plot data for a node pair N1,N2 using
	% X = squeeze(SFdata.RxTime(N1,N2,:)) - SFdata.TxTime(:,N2);
	% plot(X-mean(X),'.');		% slope of this data is frequency difference
	
	% determine round trip delay
	RTDvalid = and(FDValid,FDValid');
	RTDdata.RTDvalid(:,:,SFidx) = RTDvalid;
	RawRTD = (CurrentTT + CurrentTT').*RTDvalid;
	RTDcorrection = - ProbeDelay.*CurrentFreqDiff.*RTDvalid;
	RTD = RawRTD + RTDcorrection;
	RTDdata.RTD(:,:,SFidx) = RTD;
    
    a1 = 10.^(double(SFdata.RxSNR(:,:,SFidx)) / 20) * (double(SFdata.RxAmp(:,:,SFidx)) / 31);
    a2 = 10.^(double(SFdata.RxSNR(:,:,SFidx)) / 20) * (double(SFdata.RxAmp2(:,:,SFidx)) / 31);
    a3 = 10.^(double(SFdata.RxSNR(:,:,SFidx)) / 20) * (double(SFdata.RxAmpMax(:,:,SFidx)) / 31);

    t1 = round( (double(SFdata.RxTimeMax(:,:,SFidx)) - double(SFdata.RxTime(:,:,SFidx)))  * (32 / 5) * 10^9 );
    t2 = round( (double(SFdata.RxTimeMax(:,:,SFidx)) - double(SFdata.RxTime2(:,:,SFidx))) * (32 / 5) * 10^9 );
    
    p1 = zeros(NumNodes,NumNodes);
    p2 = zeros(NumNodes,NumNodes);
    p3 = zeros(NumNodes,NumNodes);
    
    for k1 = 1:NumNodes
        for k2 = 1:NumNodes
            if (RTDvalid(k1,k2)~=0)
                locn1 = find(prob_table1>a1(k1,k2),1,'first');
                locn2 = find(prob_table2>t1(k1,k2),1,'first');
                
                if isempty(locn1)
                    p1(k1,k2) = 0.99;
                else
                    if ~isempty(locn2)
                        p1(k1,k2) = good_prob(locn1,locn2);
                    else
                        p2(k1,k2) = good_prob(locn1,8);
                    end
                end
                
                locn1 = find(prob_table1>a2(k1,k2),1,'first');
                locn2 = find(prob_table2>t2(k1,k2),1,'first');
                
                if isempty(locn1)
                    p2(k1,k2) = 0.99;
                else
                    if ~isempty(locn2)
                        p2(k1,k2) = good_prob(locn1,locn2);
                    else
                        p2(k1,k2) = good_prob(locn1,8);
                    end
                end
                
                locn1 = find(prob_table1>a3(k1,k2),1,'first');
                
                if isempty(locn1)
                    p3(k1,k2) = 0.99;
                else
                    p3(k1,k2) = good_prob(locn1,1);
                end
            end
        end
    end
    
    RTDdata.Prob(:,:,SFidx) = p1 .* p1';

	% calculate other RTD values
	if CalcAllRTD~=0
		CurrentRxTime2 = squeeze(SFdata.RxTime2(:,:,SFidx));
		CurrentRxTimeMax = squeeze(SFdata.RxTimeMax(:,:,SFidx));
        
		CurrentTT2 = CurrentRxTime2 - TxTimeArray;
		CurrentTTMax = CurrentRxTimeMax - TxTimeArray;
		
		RTD = (CurrentTT + CurrentTT2').*RTDvalid + RTDcorrection;
		RTDdata.RTD5(:,:,1,SFidx) = RTD;
        RTDdata.Prob5(:,:,1,SFidx) = p1 .* p2';
        
		RTD = (CurrentTT + CurrentTTMax').*RTDvalid + RTDcorrection;
		RTDdata.RTD5(:,:,2,SFidx) = RTD;
        RTDdata.Prob5(:,:,2,SFidx) = p1 .* p3';
        
		RTD = (CurrentTT2 + CurrentTT2').*RTDvalid + RTDcorrection;
		RTDdata.RTD5(:,:,3,SFidx) = RTD;
        RTDdata.Prob5(:,:,3,SFidx) = p2 .* p2';
        
		RTD = (CurrentTT2 + CurrentTTMax').*RTDvalid + RTDcorrection;
		RTDdata.RTD5(:,:,4,SFidx) = RTD;
        RTDdata.Prob5(:,:,4,SFidx) = p2 .* p3';
                
		RTD = (CurrentTTMax + CurrentTTMax').*RTDvalid + RTDcorrection;
		RTDdata.RTD5(:,:,5,SFidx) = RTD;
        RTDdata.Prob5(:,:,5,SFidx) = p3 .* p3';
        
        RTD = (CurrentTT2 + CurrentTT').*RTDvalid + RTDcorrection;
        RTDdata.RTD5(:,:,6,SFidx) = RTD;
        RTDdata.Prob5(:,:,6,SFidx) = p2 .* p1';
        
        RTD = (CurrentTT2' + CurrentTTMax).*RTDvalid + RTDcorrection;
        RTDdata.RTD5(:,:,7,SFidx) = RTD;
        RTDdata.Prob5(:,:,7,SFidx) = p2' .* p3;
        
        RTD = (CurrentTT' + CurrentTTMax).*RTDvalid + RTDcorrection;
        RTDdata.RTD5(:,:,8,SFidx) = RTD;
        RTDdata.Prob5(:,:,8,SFidx) = p2' .* p3';
	end
end
