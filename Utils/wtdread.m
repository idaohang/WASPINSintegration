% Phil's wtdread

function [d, n] = wtdread(filename)
%% FUNCTION   d = wtdread(filename)  Read WASP TOA data file
%% ARGUMENT:  filename Input file name with extension .wtd
%%  RETURNS:
%%          d       Nx33 array defined as follows:
%%          column 1     sequence number
%%                 2     RX node ID
%%                 3     TX node ID
%%                 4     slot number
%%                 5     high 32-bit of RX time - in 10.24us units
%%                 6     TX time stamp - in 20.48us units
%%                 7     TOA fine time - low 16-bit of RX time in 5/32ns units
%%                 8     OLD: TOA fine time 2 
%%						 NEW: Second TOA relative to main TOA
%%                 9     OLD: TOA abs maximum
%%						 NEW: Max TOA relative to main TOA
%%                 10    OLD: spread of TOA fine
%%						 NEW: Second TOA amplitude as fraction of highest peak
%%                 11    OLD: probability of TOA fine
%%						 NEW: Max TOA amplitude as fraction of highest peak
%%                 12    OLD: spread of TOA fine 2
%%						 NEW: Main TOA amplitude as fraction of highest peak
%%                 13    RX gain
%%                 14    wideband SNR in dB
%%                 15    TOA error code
%%                 16    TX gain of RX node
%%                 17    battery voltage of RX node in mV
%%                 18    total trigger count
%%                 19    total RX error count
%%                 20    number of TX abandoned because of busy channel
%%                 21    number of TX abandoned because receive in progress
%%                 22    number of data CRC error
%%                 23    number of TX/RX slot mis-matches
%%                 24    collector's slot number
%%                 25    OLD: template index
%%					     NEW: TOA error bit 0, 1 and 2 with bit 13 ON to indicate new TOA data format
%%                 26-33 band 0-7 PSNR
%%  
%%      Phil Ho, March 2009

persistent pDir
if isempty(pDir), pDir = '.'; end
if nargin == 0
    [file, dirpath, fin] = uigetfile('*.wtd', 'Pick TOA data file', pDir, 'MultiSelect', 'on');
    if fin == 0 || isequal(dirpath, 0)
        d = [];
		n = [];
        return
    end
    pDir = dirpath;
	if ~iscell(file)
		filename = fullfile(pDir, file);
	else
		filename = cell(size(file));
		for k = 1:length(file)
			filename{k} = fullfile(pDir, file{k});
		end
	end
end

if exist('wtd') == 3
    [d, n] = wtd(filename);
    d = d';
else
    d = gtoa(filename);
	n = 0;
end

end

function d = gtoa(filename)
%% FUNCTION d = gtoa(filename)  Read TOA file
%% ARGUMENT: filename Input file name
%% WTD file consists of variable number of records in following format
%% record_header, diagnostic_info, beacon1_beacon2_..._beacon_n
%% record header is 8 byte long, where:
%%       byte 1      F4
%%       byte 2      number of beacons in record
%%       byte 3-4    RX node address
%%       byte 5-6    battery voltage in mV of RX node
%%       byte 7      TX gain in dB of RX node
%%       byte 8      collector's slot number
%% diagnostic info is 12-byte long, where
%%       byte 1-4    total number of triggers
%%       byte 5-8    total number of RX errors, bad header CRC 
%%       byte 9-10   number of TX abandoned due to busy channel
%%       byte 11-12  number of TX abandoned due to channel is receiving
%% beacon record is 24-byte long, where
%%       byte 1-16   8 16-bit word TOA results as decribed by Dave
%%              word | bits   | contents
%%                 0               TOA_fine(15:0)
%%                 1               TOA_fine(31:16)
%%                 2               TOA_coarse(15:0)
%%                 3      (3:0)    TOA_reliability : spread TOA_fine
%%                        (7:4)    TOA_reliability : Probability TOA_fine
%%                        (11:8)   TOA_reliability : spread TOA_fine2
%%                        (15:12)  wideband SNR quantized to 4 bits
%%                 4      (9:0)    TOA_fine2
%%                        (15:10)  TOA_abs_max(5:0)
%%                 5      (2:0)    TOA_abs_max(8:6)
%%                        (15:3)   Template index selected (13 bits, 0 to 8191)
%%                 6      (2:0)    subband 0 PSNR (Quantized to 3 bits)
%%                        (5:3)    subband 1 PSNR
%%                        (8:6)    subband 2 PSNR
%%                        (11:9)   subband 3 PSNR
%%                        (14:12)  subband 4 PSNR
%%                        (15)     subband 5 PSNR bit 0
%%                 7      (1:0)    subband 5 PSNR bits (2:1)
%%                        (4:2)    subband 6 PSNR
%%                        (7:5)    subband 7 PSNR
%%                        (13:8)   rx gain (bit 13 high gain/low gain, bits (12:0) gain setting
%%                        (15:14)  Error code 0 No error, 1 TOA bands in buffer is wrong, 2 TOA information/buffer timestamp mismatch, 3 no match in either buffer.
%%      byte 17-20  TX time stamp in 20.48us units
%%      byte 21-22  TX node ID
%%      byte 23     super-frame number
%%      byte 24     slot number
%%      Phil Ho, March 2009

try
    fid = fopen(filename, 'r', 'b');

    N = 100;     %% max number of rows
    M = 24;      %% number of columns
    n = 0;       %% current number of rows
    d = zeros(N, M);
    
    while ~feof(fid)
        %%
        %% Extract record header - 8 bytes
        %% uint8 is F4
        %% uint8 is number of beacons
        %% uint16 is device address of reporting node
        %% uint16 is battery voltage in mV in reporting node
        %% uint8 is current transmit gain in dB in reporting node
        %% uint8 is collector's slot number
        %%
        hdr = fread(fid, 4, 'uint16');
        if isempty(hdr), break, end
 
        if ~isequal(dec2hex(bitshift(hdr(1), -8)), 'F4')
            error('Bad TOA data file at row %d', n);
        end

        nB = bitand(hdr(1), 255);     %% number of beacons in record
        if nB > 0
            %%
            %% Read diagnostic information - 16 bytes
            %% uint32 total number of triggers
            %% uint32 total number of RX frame with bad header CRC
            %% uint16 number of TX aborted due to busy channel
            %% uint16 number of TX aborted while channel is receiving
            %% uint16 number of DSP CRC error on PHY layer
            %% uint16 number of RX error due to RX mis-matched slot numbers
            frS = fread(fid, 2, 'uint32');     %% frame status
            txS = fread(fid, 2, 'uint16');     %% TX status
            rxS = fread(fid, 2, 'uint16');     %% DSP RX error status

            %%
            %% Reading TOA from beacon - 24 bytes
            %% 8 uint16 are TOA result
            %% uint32 is TX time stamp
            %% uint16 is TX node address
            %% uint8 is sequence number
            %% uint8 is slot number
            %%
            for k = 1:nB
                rxT  = fread(fid, 8, 'uint16');
                txT  = fread(fid, 1, 'uint32');
                addr = fread(fid, 1, 'uint16');
                seq  = fread(fid, 1, 'uint8');
                slot = fread(fid, 1, 'uint8');
                n = n + 1;
                if n > N
                    d = [d; zeros(100, M)];        %% extend by 100 rows
                    N = size(d, 1);
                end
                d(n, 1)  = seq;                      %% sequence number
                d(n, 2)  = hdr(2);                   %% RX ID
                d(n, 3)  = addr;                     %% TX ID
                d(n, 4)  = slot;                     %% slot number
                d(n, 5)  = bitor(bitshift(rxT(3), 16), rxT(2));    %% RX time stamp in 10.24us units
                d(n, 6)  = txT;                      %% TX time stamp
                d(n, 7)  = rxT(1);                   %% TOA fine in 5/32ns units
                d(n, 8)  = bitand(rxT(5), 2^10-1);   %% TOA fine 2
                d(n, 9)  = bitand(bitshift(rxT(5), -10), 2^6-1);      %% TOA max
                d(n, 9)  = bitor(d(n, 9), bitshift(bitand(rxT(6), 2^3-1), 6)); 
                d(n, 10) = bitand(rxT(4), 2^4-1);    %% spread of TOA fine
                d(n, 11) = bitand(bitshift(rxT(4), -4), 2^4-1);       %% probability of TOA fine
                d(n, 12) = bitand(bitshift(rxT(4), -8), 2^4-1);       %% spread of TOA fine 2
                d(n, 13) = bitand(bitshift(rxT(8), -8), 2^5-1);       %% RX gain
                if bitget(rxT(8), 14), d(n, 13) = bitor(d(n, 13), 2^5 + 2^6); end
                d(n, 14) = bitand(bitshift(rxT(4), -12), 2^4-1)*4;    %% wideband SNR
                d(n, 15) = bitand(bitshift(rxT(8), -14), 2^2-1);      %% error code
                d(n, 16) = bitshift(hdr(4), -8);                      %% tx gain
                d(n, 17) = hdr(3);          %% battery voltage
                d(n, 18) = frS(1);          %% trigger count
                d(n, 19) = frS(2);          %% frame error count
                d(n, 20) = txS(1);          %% TX while busy count
                d(n, 21) = txS(2);          %% TX while receive count 
                d(n, 22) = rxS(1);          %% RX data CRC error count
                d(n, 23) = rxS(2);          %% TX/RX slot mis-match count
                d(n, 24) = bitand(hdr(4), 255);  %% listener's slot number
            end
        end
    end
    fclose(fid);
catch
    if fid >= 0, fclose(fid); end
    rethrow(lasterror)
end

d = d(1:n, :);

end




