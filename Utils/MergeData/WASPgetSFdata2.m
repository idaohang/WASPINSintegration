%	WASPgetSFdata2.m
%
%	Mark Hedley
%	Created: 1 April 2009
%	Modified: Jan 2010 - new TOA format, dont convert input to double
%
% This function takes data from a binary WASP file using the OFDM MAC, and the calibrated excess delay due to
% the receiver gain setting, and returns a structure organised into superframes. The first (presumably
% incomplete) and last superframes are ignored and data is kept following the first changeover in sequence number.
% The data is assumed to be in superframe order, and messages are displayed if:
%   * sequence numbers are not consecutive (warning)
%   * log data reports error (error - beacon data rejected)
%	* multiple reports of transmit timestamp, gain and slot from same transmitter not identical (error - 
%	  transmitter set to invalid)
% The function calculates the transmit time and three receive times, converts them to seconds and corrects the
% receive times for the receiver gain.


function [SFdata NodeList] = WASPgetSFdata2(LogData,ExcessDelay)

SFrollover = 256;

% correct if transposed
if size(LogData,1)==33
	LogData = LogData';
end

%% 

% Input data has following columns
InSeqNum = 1;			% sequence number
InRxID = 2;				% Rx node ID
InTxID = 3;				% Tx node ID
InSlot = 4;				% slot of transmitting node
InRxTimehigh = 5;		% high 32-bit of TOA, units of 20.48/2us
InTxTime = 6;			% transmit time stamp in units of 20.48us
InRxTimelow = 7;		% low 16-bit of TOA, units of 5/32 ns (this is a continuous binary number with InRxTimehigh
InTOA2 = 8;				% fine time for 2nd possible TOA, 15-bit signed offset to primary TOA
InTOAmax = 9;			% TOA abs maximum, 15-bit signed offset to primary TOA
InRxAmp2 = 10;			% Second TOA amplitude as fraction of highest peak
InRxAmpMax = 11;		% Max TOA amplitude as fraction of highest peak
InRxAmp1 = 12;			% First TOA amplitude as fraction of highest peak
InRxGain = 13;			% receiver gain setting (as per MAX2829 register setting)
InSNR = 14;				% wideband SNR in dB (6 bit number)
InErrCode = 15;			% TOA error code, 0 for no error
InTxGain = 16;			% TX gain of RX node
InVbat = 17;			% battery voltage of RX node in mV
InTrigCnt = 18;			% total trigger count (diagnostics)
InRxErrCnt = 19;		% total RX error count (diagnostics)
InTxAbandon = 20;		% number of TX abandoned because of busy channel (diagnostics)
InTxRxClash = 21;		% number of TX abandoned because receive in progress (diagnostics)
InCRCerr = 22;			% number of data CRC error
InSlotErr = 23;			% number of TX/RX slot mis-matches
%InCollectorSlot = 24;	% collector's slot number (i.e. slot in which node receiving this beacon normally transmits)
InErrBit = 25;			% TOA error bit 0, 1 and 2 with bit 13 ON to indicate new TOA data format
InBandPSNR = 26:33;		% band 0-7 PSNR, 2 bits each, 0 = 0-20dB, 1 = 20-30dB, 2 = 30-40dB, 3 = >40dB

% group inputs
InDiagnostics = [InTrigCnt InRxErrCnt InTxAbandon InTxRxClash InCRCerr InSlotErr];
NumBands = length(InBandPSNR);
NumDiag = length(InDiagnostics);

% calculate transmit and receive times
TxTimeIn = double(LogData(:,InTxTime))*20.48e-6;	% convert to seconds
RxTime1In = (double(LogData(:,InRxTimehigh))*double(2^16) + double(LogData(:,InRxTimelow)))*5/32*1e-9; % primary rx time in sec
RxTime2In = RxTime1In + (double(LogData(:,InTOA2)) - (2^15)*double(LogData(:,InTOA2) > (2^14)))*5/32*1e-9;
RxTimeMaxIn = RxTime1In + (double(LogData(:,InTOAmax)) - (2^15)*double(LogData(:,InTOAmax) > (2^14)))*5/32*1e-9;

% find first and last SF and eliminate as assumed partial data
SNdif = int32(LogData(2:end,InSeqNum)) - int32(LogData(1:(end-1),InSeqNum));
SNch = find(SNdif ~= 0);
if length(SNch) < 3
	error('Data too short, only %d changes in sequence number',length(SNch));
end
OutRange = (SNch(1)+1):SNch(end);
%OutLen = length(OutRange);
NumSF = length(SNch)-1;

% determine nodes used
NodeIDused = [LogData(OutRange,InRxID); LogData(OutRange,InTxID)];
NodeList = sortrows(unique(NodeIDused));
NumNodes = length(NodeList);

% generate map from Node ID to row index into NodeLocn
NodeIDmap = zeros(1,max(NodeList));
NodeIDmap(NodeList) = 1:NumNodes;

% sort data into superframes, assume received data in order
% aSeqNum = int16(zeros(NumSF,1));
% aTxValid = int16(zeros(NumSF,NumNodes));
% aTxTime = zeros(NumSF,NumNodes);
% aTxSlot = int16(zeros(NumSF,NumNodes));
% aTxGain = int16(zeros(NumSF,NumNodes));
% aVbat = int16(zeros(NumSF,NumNodes));
% aRxValid = int16(zeros(NumNodes,NumNodes,NumSF));
% aRxTime1 = zeros(NumNodes,NumNodes,NumSF);
% aRxAmp1 = int16(zeros(NumNodes,NumNodes,NumSF));
% aRxTime2 = zeros(NumNodes,NumNodes,NumSF);
% aRxAmp2 = int16(zeros(NumNodes,NumNodes,NumSF));
% aRxTimeMax = zeros(NumNodes,NumNodes,NumSF);
% aRxAmpMax = int16(zeros(NumNodes,NumNodes,NumSF));
% aRxGain = uint16(zeros(NumNodes,NumNodes,NumSF));
% aRxSNR = int16(zeros(NumNodes,NumNodes,NumSF));
% aInErrBit = int16(zeros(NumNodes,NumNodes,NumSF));
% aDiagnostics = int32(zeros(NumNodes,NumDiag,NumSF));
% aInBandPSNR = int8(zeros(NumNodes,NumNodes,NumBands,NumSF));

aSeqNum = zeros(NumSF,1,'int16');
aTxValid = zeros(NumSF,NumNodes,'int16');
aTxTime = zeros(NumSF,NumNodes);
aTxSlot = zeros(NumSF,NumNodes,'int16');
aTxGain = zeros(NumSF,NumNodes,'int16');
aVbat = zeros(NumSF,NumNodes,'int16');
aRxValid = zeros(NumNodes,NumNodes,NumSF,'int16');
aRxTime1 = zeros(NumNodes,NumNodes,NumSF);
aRxAmp1 = zeros(NumNodes,NumNodes,NumSF,'int16');
aRxTime2 = zeros(NumNodes,NumNodes,NumSF);
aRxAmp2 = zeros(NumNodes,NumNodes,NumSF,'int16');
aRxTimeMax = zeros(NumNodes,NumNodes,NumSF);
aRxAmpMax = zeros(NumNodes,NumNodes,NumSF,'int16');
aRxGain = zeros(NumNodes,NumNodes,NumSF,'uint16');
aRxSNR = zeros(NumNodes,NumNodes,NumSF,'int16');
aInErrBit = zeros(NumNodes,NumNodes,NumSF,'int16');
aDiagnostics = zeros(NumNodes,NumDiag,NumSF,'int32');
aInBandPSNR = zeros(NumNodes,NumNodes,NumBands,NumSF,'int8');
for SFidx = 1:NumSF
	
	% read data into arrays
% 	SFvalid = int16(zeros(NumNodes,NumNodes));
% 	SFslot = int16(zeros(NumNodes,NumNodes));
% 	SFtxTime = zeros(NumNodes,NumNodes);
% 	SFrxTime1 = zeros(NumNodes,NumNodes);
% 	SFrxAmp1 = int16(zeros(NumNodes,NumNodes));
% 	SFrxTime2 = zeros(NumNodes,NumNodes);
% 	SFrxAmp2 = int16(zeros(NumNodes,NumNodes));
% 	SFrxTimeMax = zeros(NumNodes,NumNodes);
% 	SFrxAmpMax = int16(zeros(NumNodes,NumNodes));
% 	SFtxGain = int16(zeros(NumNodes,NumNodes));
% 	SFrxGain = uint16(zeros(NumNodes,NumNodes));
% 	SFSNR = int16(zeros(NumNodes,NumNodes));
% 	SFVbat = int16(zeros(NumNodes,NumNodes));
% 	SFdiagnostics = int32(zeros(NumNodes,NumDiag));
% 	SFbandPSNR = int8(zeros(NumNodes,NumNodes,NumBands));
% 	SFInErrBit = int16(zeros(NumNodes,NumNodes));

	SFvalid = zeros(NumNodes,NumNodes,'int16');
	SFslot = zeros(NumNodes,NumNodes,'int16');
	SFtxTime = zeros(NumNodes,NumNodes);
	SFrxTime1 = zeros(NumNodes,NumNodes);
	SFrxAmp1 = zeros(NumNodes,NumNodes,'int16');
	SFrxTime2 = zeros(NumNodes,NumNodes);
	SFrxAmp2 = zeros(NumNodes,NumNodes,'int16');
	SFrxTimeMax = zeros(NumNodes,NumNodes);
	SFrxAmpMax = zeros(NumNodes,NumNodes,'int16');
	SFtxGain = zeros(NumNodes,NumNodes,'int16');
	SFrxGain = zeros(NumNodes,NumNodes,'uint16');
	SFSNR = zeros(NumNodes,NumNodes,'int16');
	SFVbat = zeros(NumNodes,NumNodes,'int16');
	SFdiagnostics = zeros(NumNodes,NumDiag,'int32');
	SFbandPSNR = zeros(NumNodes,NumNodes,NumBands,'int8');
	SFInErrBit = zeros(NumNodes,NumNodes,'int16');

	% check sequence number consecutive
	aSeqNum(SFidx) = LogData(SNch(SFidx)+1,InSeqNum);
	if SFidx ~= 1
		if aSeqNum(SFidx) ~= mod(aSeqNum(SFidx-1)+1,SFrollover)
			disp(sprintf('WARNING sequence number not consecutive, %d to %d', aSeqNum(SFidx-1), aSeqNum(SFidx) ));
		end
	end
	
	% put all SF data in arrays
	for Drow=(SNch(SFidx)+1):SNch(SFidx+1)
		
		% check error code
		if LogData(Drow,InErrCode)==2
% 			disp(sprintf('WARNING code 2 in logged data for SeqNum=%d, Rx=%d and Tx=%d',...
% 				OutData(Drow,OutSeqNum),OutData(Drow,OutRxID),OutData(Drow,OutTxID)));
		elseif LogData(Drow,InErrCode) ~= 0
			disp(sprintf('ERROR code %d in logged data for SeqNum=%d, Rx=%d and Tx=%d',LogData(Drow,InErrCode),...
				aSeqNum(SFidx),LogData(Drow,InRxID),LogData(Drow,InTxID)));
			continue;
		end
		if bitand(LogData(Drow,InErrBit),4) ~= 0
			disp(sprintf('No TOA in logged data for SeqNum=%d, Rx=%d and Tx=%d',...
				aSeqNum(SFidx),LogData(Drow,InRxID),LogData(Drow,InTxID)));
			continue;
		end
		
		% put data in matrices
		RxIdx = NodeIDmap(LogData(Drow,InRxID));
		TxIdx = NodeIDmap(LogData(Drow,InTxID));
		SFvalid(RxIdx,TxIdx) = 1;
		SFslot(RxIdx,TxIdx) = LogData(Drow,InSlot);
		SFtxTime(RxIdx,TxIdx) = TxTimeIn(Drow);
		SFrxTime1(RxIdx,TxIdx) = RxTime1In(Drow);
		SFrxAmp1(RxIdx,TxIdx) = LogData(Drow,InRxAmp1);
		SFrxTime2(RxIdx,TxIdx) = RxTime2In(Drow);
		SFrxAmp2(RxIdx,TxIdx) = LogData(Drow,InRxAmp2);
		SFrxTimeMax(RxIdx,TxIdx) = RxTimeMaxIn(Drow);
		SFrxAmpMax(RxIdx,TxIdx) = LogData(Drow,InRxAmpMax);
		SFtxGain(RxIdx,TxIdx) = LogData(Drow,InTxGain);
		SFrxGain(RxIdx,TxIdx) = LogData(Drow,InRxGain);
		SFSNR(RxIdx,TxIdx) = LogData(Drow,InSNR);
		SFVbat(RxIdx,TxIdx) = LogData(Drow,InVbat);
		SFdiagnostics(RxIdx,:) = LogData(Drow,InDiagnostics);
		SFbandPSNR(RxIdx,TxIdx,:) = LogData(Drow,InBandPSNR);
		SFInErrBit(RxIdx,TxIdx) = LogData(Drow,InErrBit);
	end
	
	% check consistency of valid tx times, tx gains, slot
	for iTx = 1:NumNodes
		aTxValid(SFidx,iTx) = 0;
		Vidx = find(SFvalid(:,iTx));
		if ~isempty(Vidx)
			
			% check transmit time
			if max(abs(SFtxTime(Vidx,iTx)-SFtxTime(Vidx(1),iTx))) ~= 0
				disp(sprintf('ERROR inconsistent transmit times at SeqNum=%d(SFidx=%d) for node %d',...
					aSeqNum(SFidx),SFidx,NodeList(iTx)))
				continue;
			else
				aTxTime(SFidx,iTx) = SFtxTime(Vidx(1),iTx);
			end
			
			% check transmit slot
			if max(abs(SFslot(Vidx,iTx)-SFslot(Vidx(1),iTx))) ~= 0
				disp(sprintf('ERROR inconsistent transmit slot at SeqNum=%d(SFidx=%d) for node %d',...
					aSeqNum(SFidx),SFidx,NodeList(iTx)))
				continue;
			else
				aTxSlot(SFidx,iTx) = SFslot(Vidx(1),iTx);
			end

			aTxValid(SFidx,iTx) = 1;
		end
	end

	% check consistency of Vbat - there is only one transmission of this value, so error here could only
	% result from error in matlab code that reads binary data
	for iRx = 1:NumNodes
		Vidx = find(SFvalid(iRx,:));
		if ~isempty(Vidx)
			% check transmit Vbat
			if max(abs(SFVbat(iRx,Vidx)-SFVbat(iRx,Vidx(1)))) ~= 0
				disp(sprintf('ERROR inconsistent battery voltage at SeqNum=%d(SFidx=%d) for node %d',...
					aSeqNum(SFidx),SFidx,NodeList(iRx)))
				continue;
			else
				aVbat(SFidx,iRx) = SFVbat(iRx,Vidx(1));
            end

            
            % check transmit gain
			if max(abs(SFtxGain(iRx,Vidx)-SFtxGain(iRx,Vidx(1)))) ~= 0
				disp(sprintf('ERROR inconsistent transmit gain at SeqNum=%d(SFidx=%d) for node %d',...
					aSeqNum(SFidx),SFidx,NodeList(iRx)))
				continue;
			else
				aTxGain(SFidx,iRx) = SFtxGain(iRx,Vidx(1));
			end

            
		end
	end
	
	% copy SF rx data
	aRxValid(:,:,SFidx) = SFvalid;
	aRxTime1(:,:,SFidx) = SFrxTime1;
	aRxAmp1(:,:,SFidx) = SFrxAmp1;
	aRxTime2(:,:,SFidx) = SFrxTime2;
	aRxAmp2(:,:,SFidx) = SFrxAmp2;
	aRxTimeMax(:,:,SFidx) = SFrxTimeMax;
	aRxAmpMax(:,:,SFidx) = SFrxAmpMax;
	aRxGain(:,:,SFidx) = SFrxGain;
	aRxSNR(:,:,SFidx) = SFSNR;
	aDiagnostics(:,:,SFidx) = SFdiagnostics;
	aInBandPSNR(:,:,:,SFidx) = SFbandPSNR;
	aInErrBit(:,:,SFidx) = SFInErrBit;

end

% convert flags to logical values
aRxValid = logical(aRxValid);
aTxValid = logical(aTxValid);

% correct RxTime for receiver gain
RxG = mod(aRxGain,32) + 1;
aRxTime1 = aRxTime1 - ExcessDelay(RxG);
aRxTime1 = aRxTime1.*aRxValid;
aRxTime2 = aRxTime2 - ExcessDelay(RxG);
aRxTime2 = aRxTime2.*aRxValid;
aRxTimeMax = aRxTimeMax - ExcessDelay(RxG);
aRxTimeMax = aRxTimeMax.*aRxValid;

% data returned by function
SFdata.NumSF = NumSF;
SFdata.SeqNum = aSeqNum;
SFdata.TxValid = aTxValid;
SFdata.TxTime = aTxTime;
SFdata.TxSlot = aTxSlot;
SFdata.TxGain = aTxGain;
SFdata.Vbat = aVbat;
SFdata.RxValid = aRxValid;
SFdata.RxTime = aRxTime1;
SFdata.RxAmp = aRxAmp1;
SFdata.RxTime2 = aRxTime2;
SFdata.RxAmp2 = aRxAmp2;
SFdata.RxTimeMax = aRxTimeMax;
SFdata.RxAmpMax = aRxAmpMax;
SFdata.RxGain = aRxGain;
SFdata.RxSNR = aRxSNR;
SFdata.Diagnostics = aDiagnostics;
SFdata.InBandPSNR = aInBandPSNR;
SFdata.InErrBit = aInErrBit;

