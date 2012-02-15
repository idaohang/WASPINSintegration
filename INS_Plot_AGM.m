%% User selects matlab data file to open

% cut and paste of Phil's code
pDir = '.';
[file, dirpath, fin] = uigetfile('*.mat', 'Pick an INS AGM data file', pDir);
if fin == 0 || isequal(file, 0) || isequal(dirpath, 0), return, end
pDir = dirpath;
filename = fullfile(pDir, file);

% Load the matlab file

load(filename);

%% the file contains the following data structures saved as
% save(newfile, 'INS_Gyro', 'INS_Accel', 'INS_Gtemp', ...
%              'INS_Magn', 'INS_Time', 'INS_Header', ...
%              'INS_File', 'INS_Node');

%% INS_Gyro is an n x 3 array of gyro data (x y z) in deg/sec

%% INS_Accel is an n x 3 array of accelerometer data (x y z) in milli-g
%               - multiply by 9.8/1000 to get m/s^2

%% INS_Gtemp is an n x 3 array of gyro temperature data (x y z) in C

%% INS_Magn is an n x 3 array of magnetometer data (x y z) in mGauss 
% for the 16405 and undefined units for the 1635x INS units
%               - multiply by 1/1000 to get Gauss (for 16405)

%% INS_Time is an n element array of time samples in sec

%% INS_Header
% is a structure containing the sensor's set up data
%       INS_Header.acq         is the Acquisition tag
%       INS_Header.start_ts    is the Start time stamp
%       INS_Header.gyro_offs   is the 3 value (x y z) gyro offset
%       INS_Header.accel_offs  is the 3 value (x y z) accelerometer offset
%       INS_Header.maghif      is the 3 value (x y z) magnetometer "hard
%                                 iron factor" a kind of offset
%       INS_Header.magsif      is the 3 value (x y z) magnetometer  "soft
%                iron factor" a scale on the values 2048 = 1x, 0 = 0x
%       INS_Header.adis_sr     is the ADIS internal sample rate register
%       INS_Header.adis.df     is the ADIS dynamic & filter taps register

%% INS_File       
% is a structure containing the original source file's data
%       INS_File.name       is the original file name string
%       INS_File.date       is the original file date string
%       INS_File.time       is the original file time string
%       INS_File.tim_len    is the length of the data file in seconds
%       INS_File.nodeoffs   is the initial time stamp (in seconds) of the data
%               (note the original time stamp values for sample "n" is
%                   INS_data(1,1) + INS_time.nodeoffs
%                the final time stamp would be
%                   INS_data(1,1) + INS_time.nodeoffs + INS_time.length)
%       INS_File.refoffs    is the equivalent time offset of the first sample
%               in the reference node time frame.
%               (note for the reference frame, the initial time values for 
%                sample "n" is
%                   INS_data(1,1) + INS_time.refoffs
%       INS_File.length     is the number of data samples in the data array
%       INS_File.comment    is a comment string to preserve critical data about the matlab preprocessing
%           typically it is 'OK'; but some data arrays had recording
%           problems and the first data samples are spurious. In this case
%           the string might be 'Trimmed 100 samples'

%% INS_Node
% is a structure containing information about the source nodes sensor
%       INS_Node.nodeno is the sensor's node number e.g. 218 or 224
%       INS_Node.sensno identifies the sensor daughter board number (1 to 5)
%       INS_Node.sensdev identifies the ADIS sensor e.g. 16350, 16355 etc.

%% Plot the INS data

% Get number of samples
nI = INS_File.length;
% Get time and delta time
t = INS_Time;

dt = (t(nI)-t(1))/nI;

figure(1);
hold off;
plot(t,INS_Gyro(:,1),'r');
hold on;
plot(t,INS_Gyro(:,2),'g');
plot(t,INS_Gyro(:,3),'b');
ylabel('X, Y, Z (deg/s)');
pName = horzcat('Gyro X, Y, Z vs Time - ',file);
title(pName);
legend('X', 'Y', 'Z');

figure(2);
hold off;
plot(t,INS_Accel(:,1),'r');
hold on;
plot(t,INS_Accel(:,2),'g');
plot(t,INS_Accel(:,3),'b');
ylabel('X, Y, Z (mG)');
pName = horzcat('Accelerometer X, Y, Z vs Time - ',file);
title(pName);
legend('X', 'Y', 'Z');

figure(3);
hold off;
plot(t,INS_Magn(:,1),'r');
hold on;
plot(t,INS_Magn(:,2),'g');
plot(t,INS_Magn(:,3),'b');
ylabel('X, Y, Z (mGauss)');
pName = horzcat('Magnetometer X, Y, Z vs Time - ',file);
title(pName);
legend('X', 'Y', 'Z');
