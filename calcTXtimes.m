

% Script puts the local transmission times into the excel worksheets with
% WASP position

clear all;
files = 15;

offset220 = zeros(files,1);
offset315 = zeros(files,1);
offset315to220 = zeros(files,1);

for i=1:files
    
   load(['WASP_INS_Data/20111220/WASP20Dec2011_rd_' num2str(i) '.mat']);
   offset220(i) = mode(TxTimesLocl(:,9)-TxTimesLocl(:,5));
   offset315(i) = mode(TxTimesLocl(:,9)-TxTimesLocl(:,10));
   offset315to220(i) = mode(TxTimesLocl(:,10)-TxTimesLocl(:,5));
   
   positions = csvread(['WASP_INS_Data/20111220/PosWASP21Dec' num2str(i) '.csv']); 
   
   % Check for missing data (based on time == 0) and ensure it's nan instead of zero
   positions(positions(:,4) == 0, 1:8) = nan;
   positions(positions(:,11) == 0, 9:15) = nan;
   
   % Make the adjustment to the time
   positions(:,4) = positions(:,4) - offset220(i);
   positions(:,11) = positions(:,11) - offset315(i);  
   
   csvwrite(['WASP_INS_Data/20111220/PosWASP21Dec' num2str(i) 'A.csv'], positions);
   
   clear positions; 
end