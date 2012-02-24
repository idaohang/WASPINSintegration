% Process raw data files and generate ranges and save them for post
% processing.
%
% Various file names are passed to the TrialMerge.m function,
% which are self explanatory. This function can handle multiple raw data
% files and removes any redundant data and generates unique range
% measurements.
%
% NodesToRemove variabile consists of node IDs of the nodes
% whose measuements are not needed.

function MergeData()

clc

addpath(genpath('Utils'));
BaseDataDir = '';

DIM = 2;                                  % two or three dimensional localization
REQUIRE_MAC_ID_TO_NODE_ID_CONVERSION = 1; % 1 - need conversion from MAC_ID to Node_ID, 0 - otherwise
AIS_OR_ADHOC = 0;                         % 1 - data in AIS mode, hence expects TDMA file

% % Marsfield Jan 22, 2012
CollectionNodes = 1;
SelectDataSet = 1;
SurveyFileName = [BaseDataDir 'WASP_INS_data/20120222/survey.csv'];
NodeDelayFileName = 'NodeDelay.csv';
FileNames = {[BaseDataDir 'WASP_INS_data/20120222/WASP22Feb' num2str(SelectDataSet) '.wtd']};
RangeDataFileName = [BaseDataDir 'WASP_INS_data/20120222/WASP20Dec2011_rd_' num2str(SelectDataSet)];
NodesToRemove = [];


% % % Marsfield Dec 20, 2011
% CollectionNodes = 1;
% SelectDataSet = 14;
% SurveyFileName = [BaseDataDir 'Marsfield_2011_12_20/Survey_20Dec2011.csv'];
% NodeDelayFileName = 'NodeDelay.csv';
% FileNames = {[BaseDataDir 'Marsfield_2011_12_20/WASP20Dec' num2str(SelectDataSet) '.wtd']};
% RangeDataFileName = [BaseDataDir 'Marsfield_2011_12_20/WASP20Dec2011_rd_' num2str(SelectDataSet)];
% NodesToRemove = [];


% % % % NSW FB data Nov 15, 2010
% CollectionNodes = 1;
% SelectDataSet = 6;
% SurveyFileName = [BaseDataDir 'NSWFB_15Nov2010/NSWFB_15Nov2010_Survey.csv'];
% NodeDelayFileName = 'NodeDelay.csv';
% FileNames = {[BaseDataDir 'NSWFB_15Nov2010/WASP15Nov' num2str(SelectDataSet) '.wtd']};
% RangeDataFileName = [BaseDataDir 'NSWFB_15Nov2010/WASP15Nov2010_rd_' num2str(SelectDataSet)];
% NodesToRemove = [];

% % % Marsfield Feb 08, 2011
% CollectionNodes = {'256', '301', '323'};
% SelectDataSet = 1;
% SurveyFileName = [BaseDataDir 'Marsfield_2011_02_08/Sathyan_8Feb2011.csv'];
% NodeDelayFileName = 'NodeDelay.csv';
% %FileNames = {[BaseDataDir 'AIS_Netball_2010_10_20/Raw20Oct1.wtd']};
% RangeDataFileName = [BaseDataDir 'Marsfield_2011_02_08/WASP08Feb2011_rd_' num2str(SelectDataSet)];
% NodesToRemove = [];

% % AIS Netball Oct 20, 2010
% CollectionNodes = 1;
% SelectDataSet = 1;
% SurveyFileName = [BaseDataDir 'AIS_Netball_2010_10_20/TDMA.csv'];
% NodeDelayFileName = [BaseDataDir 'AIS_Netball_2010_10_20/NodeDelayCal.csv'];
% FileNames = {[BaseDataDir 'AIS_Netball_2010_10_20/Raw20Oct1.wtd']};
% RangeDataFileName = [BaseDataDir 'AIS_Netball_2010_10_20/WASP20Oct2010_rd_' num2str(SelectDataSet)];
% NodesToRemove = [];

% % AIS Canberra 23 August 2011 
% CollectionNodes = 1;
% SelectDataSet = 1;
% RangeDataFileName = [BaseDataDir 'AIS_Canberra_2011_08_23/WASP23Aug2011_rd_' num2str(SelectDataSet)];
% SurveyFileName = [BaseDataDir 'AIS_Canberra_2011_08_23/TDMA.csv'];
% NodeDelayFileName = [BaseDataDir 'AIS_Canberra_2011_08_23/NodeDelayX.csv'];
% FileNames = {[BaseDataDir 'AIS_Canberra_2011_08_23/WASP23Aug' num2str(SelectDataSet) '.wtd']};
% NodesToRemove = [];

% % AIS Canberra 24 August 2011 
% CollectionNodes = 1;
% SelectDataSet = 1;
% SurveyFileName = [BaseDataDir 'AIS_Canberra_2011_08_24/TDMA.csv'];
% NodeDelayFileName = [BaseDataDir 'AIS_Canberra_2011_08_24/NodeDelayX.csv'];
% FileNames = {[BaseDataDir 'AIS_Canberra_2011_08_24/WASP24Aug' num2str(SelectDataSet) '.wtd']};
% RangeDataFileName = [BaseDataDir 'AIS_Canberra_2011_08_24/WASP24Aug2011_rd_' num2str(SelectDataSet)];
% NodesToRemove = [];


% %% AIS Canberra 24 August 2011 
% CollectionNodes = 1;
% SelectDataSet = 3;
% SurveyFileName = [BaseDataDir 'AIS_Canberra_2011_11_25/AIS_Survey_25Nov2011.csv'];
% NodeDelayFileName = 'NodeDelay.csv';
% FileNames = {[BaseDataDir 'AIS_Canberra_2011_11_25/w307-20111125-0000' num2str(SelectDataSet) '.wtd']};
% RangeDataFileName = [BaseDataDir 'AIS_Canberra_2011_11_25/WASP25Nov2011_rd_' num2str(SelectDataSet) '_6anc'];
% NodesToRemove = [400, 407, 408, 427, 250, 261, 233, 243, 248, 228];

%% Number of files should be the same as the number of collection nodes
NumFiles = length(CollectionNodes);
if (NumFiles > 1)    
    FileNames = cell(NumFiles,1);   
    % Generate filenames according to the selected data set
    for k1 = 1:NumFiles
        FileNames{k1} = [BaseDataDir 'Marsfield_2011_02_08/w' CollectionNodes{k1} '-20110208-0000' num2str(SelectDataSet) '.wtd'];
    end
end

%% The name of the RxGaindelay file should go here - file with additional
%% receiver delay due to gain
FileNameRxGainDelay = 'RxGainDelayFile.csv';

TrialMerge(FileNameRxGainDelay,NodeDelayFileName,SurveyFileName,FileNames,...
    RangeDataFileName,NodesToRemove,DIM,REQUIRE_MAC_ID_TO_NODE_ID_CONVERSION,AIS_OR_ADHOC);

