
%function [ r, Cn2b, timestamps_IMU, WASPloc, timestamps_WASP  ] = main( INSfile, loss_freq, loss_period )
function [  ] = main( INSfile )

% main function that runs an INS / WASP Kalman filter integration
% Do a 'clear all' first or sometimes you will get an error
clearvars -except INSfile loss_freq loss_period 
global a v r rpy Cn2b ba bg dr dv fs fsc ws wsc epsilon WASPloc timestamps_IMU...
 timestamps_WASP nS nS_filter front Sk innov ms mn_cal gn_cal Pk msc bm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   LOADING DATA

% Read imu data stored as Nx13 array of TsGxGyGzAxAyAzTxTyTzMxMyMz
if nargin < 1
    [imu,imu_p,filename,dirpath] = wsdread();
else
    [imu,imu_p,filename,dirpath] = wsdread(INSfile);
end;

% Map the  INS wsd file to the correct WASP file
mapfile = fopen(fullfile(dirpath,'filemapping.csv')); 
filemapping = textscan(mapfile,'%s %s %s %d %f %f %f %f %f %f', 'delimiter',',');
fclose(mapfile); 
index = find(ismember(filemapping{1}, filename)==1);
WASPdatafile = strcat(dirpath,filemapping{2}(index));
col = filemapping{4}(index); % Which column to use out of the posWASP csv files
front = [filemapping{5}(index); filemapping{6}(index); filemapping{7}(index)];
ms_bias = [filemapping{8}(index); filemapping{9}(index); filemapping{10}(index)];


% Load the range data
WASPrangefile =  fullfile(dirpath,char(filemapping{3}(index)));
R = load(WASPrangefile);
node = strtok(filename,'-');    
node = str2double(node(2:end));
ranges = (squeeze(R.RangeData(:,node==R.NodeList,:)) + squeeze(R.RangeData(node==R.NodeList,:,:))) ./2 ;
ranges(R.MobileNodeLocn,:) = []; %remove ranges to mobile nodes
timestamps_ranges = R.TxTimesLocl(:,node==R.NodeList);

% Load the survey file
survey = sortrows( csvread(fullfile(dirpath,'survey.csv')), 1);
survey(:,1) = [];
% Ensure surveyed positions are in same order as ranges
clear R;

% Load the WASP position data in format "TsXposYpos"
WASPdata = csvread(WASPdatafile{1},0,col);
WASPdata(:,col+1:size(WASPdata,2)) = []; % drop fields not needed
WASPdata(1:4,:) = []; % drop first 4 obs
WASPdata( isnan(WASPdata(:,1)) , :) = []; % drop NaN times
WASPdata(:,2) = -WASPdata(:,2); % flip the x-axis on WASP to go from survey to map coordinates

% Synchronise the start of IMU and WASP
imu(1:find(imu(:,1) < WASPdata(1,1),1,'last'),:) = [];
WASPdata(1:find(WASPdata(:,1) < imu(1,1),1,'last'), :) = [];

% Synchronise the end of IMU and WASP
imu(find(imu(:,1) > WASPdata(end,1),1,'first'):end,:) = [];
WASPdata(find(WASPdata(:,1) > imu(end,1),1,'first'):end,:) = [];

% Create some more meaningful names, accel for integration and orientation
% are separated with different high freq denoising constants
ws = imu(:,2:4)';   % Gyro in s-frame
fs = -imu(:,5:7)';   % Accel in s-frame (for integration)
ms = imu(:,11:13)'; % Mag in s-frame
WASPloc = WASPdata(:,2:3)';
fs_var = LocalVar(82, (imu(:,5).^2 + imu(:,6).^2 + imu(:,5).^2).^0.5);
% Local variange of the accelerometer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SETUP

samples = 820; % How many samples used to estimate initial orientation, 
% (about 820 samples in 1 second)

% Timestamps
timestamps_IMU = imu(:,1)';% - imu(1,1);
timestamps_WASP = WASPdata(:,1)';% - WASPdata(1,1);
nS = size(timestamps_IMU,2);
nS_filter = size(timestamps_WASP,2);
dt = (timestamps_IMU(end)- timestamps_IMU(1))/(nS-1);
dt_WASP = (timestamps_WASP(end)- timestamps_WASP(1))/(nS_filter-1);

% Allocate error feedback and INS output variables
[a, v, r, rpy] = deal(zeros(3,nS+1)); % INS output
[Cn2b] = deal(zeros(3,3,nS+1)); % INS orientation matrix
[msc, wsc, fsc] = deal(zeros(3,nS)); % Bias corrected gyro and accel in s-frame
[bm, ba, bg, epsilon, dv, dr] = deal(zeros(3,nS+1)); % filter error estimates
[innov, Sk] = deal(zeros(9,nS+1));       % filter residuals + std dev
[Pk] = deal(zeros(15,nS+1));             % filter std dev

% Initialise INS position using first WASP observation 
r(:,1) = [WASPloc(:,1);0]; % Assume z = 0;

% Environment calibrated magnetic vector and assumed gravity vector in
% n-frame (refer to plotmap function for the n-frame orientation)
mn_cal = RPY2DCM(deg2rad([0; 0; 207.5]))*[242; -5.3; 515];
% This is calculated from http://www.ngdc.noaa.gov/geomagmodels/struts/calcIGRFWMM
% then uses http://googlecompass.com/ to estimate the rotation of the
% n-frame from north
gn_cal = [0; 0; 1000]; % Doesn't seem to matter which one is used

% Remove previously estimated magnetometer bias
for t=1:nS
    ms(:,t) = ms(:,t) - ms_bias;
end;

% Estimate initial orientation from observed magnetic vector and gravity 
gn = mean(fs(:,1:samples),2);
mn = mean(ms(:,1:samples),2);
mn_cal_perp = mn_cal - dot(mn_cal,gn_cal)/dot(gn_cal,gn_cal)*gn_cal;
mn_perp = mn - dot(mn,gn)/dot(gn,gn)*gn; % Component of magnetic vector in xy plane
Cb2n = wahba([gn_cal,mn_cal_perp],[gn,mn_perp]); 
Cn2b(:,:,1) = Cb2n';
Cs2b = eye(3); % Assume body and sensor framed are aligned


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESSING LOOP

tf = 1;     % timestep for the filter
last_t = 1; % for averaging accelerations over a WASP filter update cycle
for t=1:nS  % timestep for the INS
    
    %bg(:,tf) = [0.42; 0.23; -0.1]; % manually tuned bias correction, for testing 315 circle
    
    % Apply bias corrections   
    msc(:,t) = ms(:,t) - bm(:,t);    % bias corrected magnetometer in s-frame
    wsc(:,t) = ws(:,t) - bg(:,t);    % bias corrected gyro in s-frame
    fsc(:,t) = fs(:,t) - ba(:,t);    % bias corrected accel in s-frame
    wb = Cs2b*wsc(:,t); % bias corrected gyro in b-frame
    fb = Cs2b*fsc(:,t); % bias corrected accel in b-frame
  
    % INS integration
    [a(:,t+1), v(:,t+1), r(:,t+1), Cb2n]...
        = INS(v(:,t), r(:,t), Cb2n, dv(:,t), dr(:,t), epsilon(:,t), wb, fb, dt, gn_cal); 
    Cn2b(:,:,t+1) = Cb2n';
    rpy(:,t+1) = rad2deg(DCM2RPY(Cb2n')); 
    
    % Error estimation using KF
    if( tf <= numel(timestamps_WASP) && t < numel(timestamps_IMU) && timestamps_WASP(tf) < timestamps_IMU(t+1) )  
        % WASP available
        
        dr_wasp = r(:,t+1) - [WASPloc(:,tf); 0]; % augment by assuming vertical movement is wrong
        if (tf > 1), dt_WASP = timestamps_WASP(tf)- timestamps_WASP(tf-1); end;
        fs_mean = mean(fsc(:,last_t:t),2); % Could do something more sophisticated here, ie. integration
        
        [ ba_error, bg_error, dv(:,t+1), dr(:,t+1), epsilon(:,t+1), innov(1:3,t+1), Sk(1:3,t+1), Pk(:,t+1) ] ...
            = WASPfilter( Cs2b, Cn2b(:,:,last_t)', fs_mean, dt_WASP, dr_wasp );
        
        ba(:,t+1) = ba(:,t) + ba_error;
        bg(:,t+1) = bg(:,t) + bg_error;
        bm(:,t+1) = bm(:,t);
                
        tf = tf + 1;        % Increment filter time step
        t = t-1;
        last_t = t;
        
    else % No WASP update available
              
        [ bm_error, bg_error, epsilon(:,t+1), innov(4:9,t+1), Sk(4:9,t+1), Pk(7:12,t+1) ] ...
            = IMUfilter( Cs2b, Cb2n, fsc(:,t), msc(:,t), dt, mn_cal, gn_cal, fs_var(t) );
        
        ba(:,t+1) = ba(:,t);
        bm(:,t+1) = bm(:,t) + bm_error;
        bg(:,t+1) = bg(:,t) + bg_error;
   
   end;
    
end;  


