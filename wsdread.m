function [c, p, file, dirpath] = wsdread(filename, arg)
%% FUNCTION   [c, p] = wsdread(filename)  Read WASP INS data file
%% ARGUMENTS: filename Input file name with extension .wsd
%%            arg      Extra argument
%%  RETURNS:
%%          c       Nx13 array of TsGxGyGzAxAyAzTxTyTzMxMyMz
%%          p       Embedded parameters as a 10x1 or 16x1 array, where
%%          p(1)    Acquisition tag
%%          p(2)    Start time stamp
%%          p(3)    X gyro offset
%%          p(4)    Y gyro offset
%%          p(5)    Z gyro offset
%%          p(6)    X accel offset
%%          p(7)    Y accel offset
%%          p(8)    Z accel offset
%% For ADIS-16350
%%          p(9)    ADIS internal sample rate register
%%          p(10)   ADIS dynamic & filter taps register
%% For ADIS-16405
%%          p(9)    x magn hard-iron factor
%%          p(10)   y magn hard-iron factor
%%          p(11)   z magn hard-iron factor
%%          p(12)   x magn soft-iron factor
%%          p(13)   y magn soft-iron factor
%%          p(14)   z magn soft-iron factor
%%          p(15)   ADIS internal sample rate register
%%          p(16)   ADIS dynamic & filter taps register
%%
%%      Phil Ho June 2009

persistent pDir bHack dSpike

if isempty(pDir), pDir = '.'; end
if isempty(bHack), bHack = 0; end
if isempty(dSpike), dSpike = 0; end

if nargin >= 1 && ischar(filename)
    if isequal(filename, 'TsHack')
        if nargin == 2 && isnumeric(arg), bHack = arg; end
        c = bHack;
        p = dSpike;
        return
    elseif isequal(filename, 'Despike')
        if nargin == 2 && isnumeric(arg), dSpike = arg; end
        c = dSpike;
        p = bHack;
        return
    else 
         % Added so it will work when given a file name
        parts = textscan(filename,'%s','delimiter','\\');
        parts = textscan(char(parts{1}(end)),'%s','delimiter','.');
        file = char(parts{1}(1));
    end
end

%%
show = 0;                               %% to plot the data ?

    
if nargin == 0, show = 1; end
if nargin == 1 && isnumeric(filename)
    show = filename;
    filename = [];
    
end

if nargin == 0 || isempty(filename)
    [file, dirpath, fin] = uigetfile('*.wsd', 'Pick an INS data file', pDir);
    if fin == 0 || isequal(file, 0) || isequal(dirpath, 0)
        c = [];
        p = [];
        return
    end
    pDir = dirpath;
    filename = fullfile(pDir, file);
    %show = 1;
end

[d, r] = gins(filename, bHack);

%%
%% Process returned variables
%%
if nargout == 0 || show
   %% View G-A-T-M INS data
   if dSpike && exist('despike', 'file') == 2, d(:, 2:end) = despike(d(:, 2:end)); end
   %%
   %% ADIS's time base of 0.61035ms, internal sampling period
   %% is (Ns + 1) * timebase, where NS is the content of r(8).
   %%
   if bHack
      t = (0:size(d, 1)-1) * (bitand(r(9), 127) + 1) * 0.61035e-3;
   else
      t = d(:, 1) - d(1, 1);           %% time axis from FPGA time stamp
   end
   [v, file] = fileparts(filename);

   subplot(4, 1, 1), plot(t, d(:, 5:7)), grid
   ylabel('Accel');
   title(sprintf('%s: T%d', file, r(1)));
   subplot(4, 1, 2), plot(t, d(:, 2:4)), grid
   ylabel('Gyro');
   subplot(4, 1, 3), plot(t, d(:, 11:13)/256), grid
   ylabel('Magneto');
   subplot(4, 1, 4), plot(t, d(:, 8:10)), grid
   ylabel('Temp');
   xlabel('time in secs');
   subplot
end
   
if nargout >= 1, c = d; end     %% return data only
if nargout >= 2, p = r; end     %% return parameters
     
end

function [c, p] = gins(filename, hack)
%% FUNCTION [c, p] = gins(filename)  Read INS file
%% ARGUMENT: filename Input file name
%%           hack     if 1 try to work around time stamp bug 
%%  RETURNS:
%%          c       Nx13 array of TsGxGyGzAxAyAzTxTyTzMxMyMz
%%          p       Embedded parameters as a 10x1 or 16x1 array, where
%%          p(1)    Acquisition tag
%%          p(2)    Start time stamp
%%          p(3)    X gyro offset
%%          p(4)    Y gyro offset
%%          p(5)    Z gyro offset
%%          p(6)    X accel offset
%%          p(7)    Y accel offset
%%          p(8)    Z accel offset
%% For ADIS-16350
%%          p(9)    ADIS internal sample rate register
%%          p(10)   ADIS dynamic & filter taps register
%% For ADIS-16405
%%          p(9)    x magn hard-iron factor
%%          p(10)   y magn hard-iron factor
%%          p(11)   z magn hard-iron factor
%%          p(12)   x magn soft-iron factor
%%          p(13)   y magn soft-iron factor
%%          p(14)   z magn soft-iron factor
%%          p(15)   ADIS internal sample rate register
%%          p(16)   ADIS dynamic & filter taps register
%%
%%      Phil Ho September 2008
    
try
    fid = fopen(filename, 'r', 'b');

    %%
    %% Extract magic signature & data type
    %%
    t = fread(fid, 2, 'uint16');
    m = dec2hex(bitshift(t(1), -8));
    if ~(isequal(m, 'F4') || isequal(m, 'F5')), error 'Bad INS data file', end
    om = isequal(m, 'F4');         %% ADIS 16350 INS module
    if om
        aScale = 2.522;            %% mg/LSB for ADIS16355
        gScale = 0.07326;          %% deg/s per LSB with 300deg range
        tScale = 0.1453;           %% deg per LSB with 25deg offset
        mScale = 1;                %% no scale for magnetometer
    else
        aScale = 10/3;             %% mg/LSB for ADIS16405
        gScale = 0.05;             %% deg/s per LSB with 300deg range
        tScale = 0.14;             %% deg/LSB with 25deg offset
        mScale = 0.5;              %% 0.5 mgauss/LSB
    end

    type = bitand(t(1), 255);
    if type == 5           %% INS data
        if t(2) ~= 32 && t(2) ~= 44, error 'Bad INS data file', end

        %%
        %% Read frame size and tag as 32-bit unsigned interger
        %%
        s1 = fread(fid, 3, 'uint32');
        if s1(1) ~= 32 && s1(1) ~= 24, error 'Bad INS FIFO frame'; end
        dSize = s1(1)/2;            %% number of words per FIFO data frame
        %%
        %% And the ADIS parameters as 16-bit values
        %%
        if om
            s2 = fread(fid, 8, 'int16');
        else
            s2 = fread(fid, 14, 'int16');
        end
        s = [s1(2:end); s2];
        rg = bitshift(s(10), -8);
        if rg == 1
            gScale = gScale/4;      %% gyro scale for 75deg/sec
        elseif rg == 2
            gScale = gScale/2;      %% gyro scale for 150deg/sec
        end
        
        %%
        %% work out number of samples in file
        %%
        pos = ftell(fid);
        fseek(fid, 0, 'eof');
        nS = (ftell(fid) - pos)/dSize/2;    %% number of samples
        fseek(fid, pos, 'bof');
        
        %%
        %% Reading sample data
        %%
        d = zeros(nS, 13);
        
        if om
            %%
            %% ADIS-16350
            %%
            for k = 1:nS
                r = fread(fid, dSize, 'uint16');
                v = unpack_350(r, gScale, aScale, tScale);
                d(k, :) = v';
            end
        else
            %%
            %% ADIS-16405
            %%
            for k = 1:nS
                r = fread(fid, dSize, 'uint16');
                v = unpack_405(r, gScale, aScale, tScale, mScale);
                d(k, :) = v';
            end
        end
        fclose(fid);
        
        %%
        %% XXXX HACK - a step change in FPGA time stamp at start
        %% that does not correspond to the start time stamp in s(2)
        %%
        if hack
            t0 = s(2)*20.48e-6;     %% start time stamp
            for k = 1:size(d, 1)
                if d(k, 1) - t0 > 0 && d(k, 1) - t0 < 0.001, break; end
            end
            if k > 1 && k < size(d, 1), d = d(k:end, :); end
        end
        
        if nargout == 1
            c = d;
        elseif nargout > 1
            c = d;
            p = s;
        end
    else
        error 'Invalid INS data'
    end
catch
    if fid >= 0, fclose(fid); end
    rethrow(lasterror)
end

end

function y = unpack_350(x, gs, as, ts)
%% FUNCTION y = unpack_350(x, gs, as)  Convert to INS sample data from FPGA
%% packed data frame.
%% ARGUMENTS:   x   Raw sample 16x1 array of uint16
%%              gs  Gyro scale
%%              as  Accel scale
%%              ts  Temperature scale
%% RETURNS:     INS samples 13x1 array of TsGxGyGzAxAyAzTxTyTzMxMyMz

%%
%% Decode INS ADIS-16350 data frame - 16 words - from FPGA packing
%%
y = zeros(13, 1);
y(1) = (x(2)*2^16 + x(1)) * 20.48e-6;               %% time stamp
y(2) = sxt(x(3), 14) * gs;                          %% gs
y(3) = sxt(x(4), 14) * gs;                          %% gy
y(4) = sxt(x(5), 14) * gs;                          %% gz
y(5) = sxt(x(6), 14) * as;                          %% ax
y(6) = sxt(x(7), 14) * as;                          %% ay
y(7) = sxt(x(8), 14) * as;                          %% az
y(8) = 25 + sxt(x(9), 12) * ts;                     %% tx                        
y(9) = 25 + sxt(x(10), 12) * ts;                    %% ty
y(10) = 25 + sxt(x(11), 12) * ts;                   %% tz
t1 = uint32(x(15)); t2 = uint32(x(16));             %% mx
y(11) = sxt(bitor(t1, bitshift(bitand(t2, 255), 16)), 24);   
t1 = uint32(x(12)); t2 = uint32(x(13));             %% my
y(12) = sxt(bitor(t1, bitshift(bitand(t2, 255), 16)), 24);
t1 = uint32(x(14));                                 %% mz
y(13) = sxt(bitor(bitshift(t1, 8), bitand(bitshift(t2, -8), 255)), 24);

end

function y = unpack_405(x, gs, as, ts, ms)
%% FUNCTION y = unpack_405(x, gs, as, ms)  Convert to INS sample data from FPGA
%% packed data frame.
%% ARGUMENTS:   x   Raw sample 12x1 array of uint16
%%              gs  Gyro scale
%%              as  Accel scale
%%              ts  Temperature scale
%%              ms  Magnetometer scale
%% RETURNS:     INS samples 13x1 array of TsGxGyGzAxAyAzTxTyTzMxMyMz

%%
%% Decode INS ADIS-16405 data frame - 12 words - from FPGA packing
%%
y = zeros(13, 1);
x(3:end) = bitand(x(3:end), 16383);                 %% mask out bits 15-16
y(1) = (x(2)*2^16 + x(1)) * 20.48e-6;               %% time stamp
y(2) = sxt(x(3), 14) * gs;                          %% gs
y(3) = sxt(x(4), 14) * gs;                          %% gy
y(4) = sxt(x(5), 14) * gs;                          %% gz
y(5) = sxt(x(6), 14) * as;                          %% ax
y(6) = sxt(x(7), 14) * as;                          %% ay
y(7) = sxt(x(8), 14) * as;                          %% az
y(8) = 25 + sxt(x(12), 12) * ts;                    %% tx                        
y(9) = y(8);                                        %% ty
y(10) = y(8);                                       %% tz
y(11) = sxt(x(9), 14) * ms;                         %% mx   
y(12) = sxt(x(10), 14) * ms;                        %% my
y(13) = sxt(x(11), 14) * ms;                        %% mz

end

function y = sxt(x, n)
%% FUNCTION y = sxt(x, n) Sign extension of an integer
%% ARGUMENTS:   x Input integer
%%              n Number of bits
%% RETURNS      Signed extended number

if x >= 2^(n - 1)
    y = x - 2^n;
else
    y = x;
end

end


