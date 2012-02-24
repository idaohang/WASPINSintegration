% Plots the range error CDF.
%
% Input:
%   RangeData - Consists of range between various pairs of nodes
%   RSSvalues - Corresponding RSS values.
%   rxSNR     - Corresponding received SNR.
%   dist1     - Distance between pairs of nodes
%
% Output:
%   rerr_ret  - Range error calculated for each measured ranges.
%
%

function rerr_ret = plot_range_consistency(RangeData, RSSvalues, rxSNR, dist1)

if nargin == 0
    RangeDataFileName = 'data/WASP05Mar23_range_data_3';
    load(RangeDataFileName);
end

dist = dist1;

MINRSS = -85;
MINSNR = 36;

numSF = size(RSSvalues,3);
numNodes = size(RSSvalues,2);

% range error distribution of individual nodes
RERRTHRESHOLD = 10;
largeRErr = zeros(numNodes,numNodes);
for Nidx = 1:numNodes
    rd = RangeData(:,Nidx,:);    
    RSSValid = (RSSvalues(:,Nidx,:) >= MINRSS) & (RSSvalues(:,Nidx,:) ~= 0);    
    rd = rd .* RSSValid;    
    SNRValid = rxSNR(:,Nidx,:) >= MINSNR;    
    rd = rd .* SNRValid;    
    rd = squeeze(rd);    
    for Nidx1 = 1:numNodes
        rd1 = rd(Nidx1,:);        
        if all(rd1==0)
            continue;
        end
        rd1(rd1==0) = [];
        rErr = abs(rd1 - dist(Nidx,Nidx1));        
        largeRErr(Nidx,Nidx1) = sum(rErr > RERRTHRESHOLD);
    end
end

% collect all the pairs of ranges in a superframe
RangeDataAll = zeros(numel(dist), numSF);
for SFidx = 1:numSF
    rd = RangeData(:,:,SFidx);    
    RSSValid = (RSSvalues(:,:,SFidx) >= MINRSS) & (RSSvalues(:,:,SFidx) ~= 0);    
    rd = rd .* RSSValid;    
    SNRValid = rxSNR(:,:,SFidx) >= MINSNR;    
    rd = rd .* SNRValid;    
    RangeDataAll(:,SFidx) = rd(:);
end
RangeDataAll(RangeDataAll==0) = Inf; % set all the zero values to Inf

% % plot the range vs true range for all the node pairs, note that this
% % plots node i - node j pair for i = j as well
% figure(10)
% clf
% hold all
% for SFidx = 1:numSF
%     plot(1:numel(dist), RangeDataAll(:,SFidx), '.');
% end
% dist(dist==0) = Inf;
% plot(1:numel(dist), dist(:), 's')

% find the range measurement error
dist = dist(:);
rerr = zeros(numel(dist),numSF);

cnt = 0;
for SFidx = 1:numSF
    ranges = RangeDataAll(:,SFidx);    
    if all(ranges==Inf)
        continue;
    else
        rd = ranges - dist;
        %rd(rd==0) = 0;
        
        % remove the NaN values due to node i - node j (i = j) and replace
        % it with Inf        
        nanIdx = isnan(rd);
        rd(nanIdx) = Inf;
        cnt = cnt + 1;
        rerr(:,cnt) = rd;
    end
end
rerr = rerr(:,1:cnt);

% remove all the Inf values and find max and min error and plot range error
% distribution
rerr1 = rerr(:);
rerr1(abs(rerr1)==Inf) = [];
maxerr = max(rerr1);
minerr = min(rerr1);
errrange = minerr:(maxerr-minerr)/100:maxerr;
N = histc(rerr1, errrange);
N = N/sum(N);

% % Filtered range error pdf - for patent, adapted for a particular data set
% windowsize = 10;
% n1 = filter(ones(1,windowsize)/windowsize,1,N);
% figure(30)
% clf
% plot(errrange,n1)
% hold on
% xmin = -0.5; xmax = 3; ymin = 0.0001; ymax = 0.05; npts = 100;
% plot(xmin*ones(1,npts),linspace(ymin,ymax,npts),'r')
% plot(linspace(xmin,xmax,npts),ymax*ones(1,npts),'r')
% plot(xmax*ones(1,npts),linspace(ymin,ymax,npts),'r')
% plot(linspace(minerr,xmin,npts),ymin*ones(1,npts),'r')
% plot(linspace(xmax,maxerr,npts),ymin*ones(1,npts),'r')
% plot(linspace(-0.975,0.97,100), linspace(0.0001,0.07,100),'k')
% plot(linspace(0.97,3.625,100), linspace(0.07,0.0001,100),'k')
% 
% axis([-2 8 0 0.08])
% xlabel('Range error (m)','FontSize',12)
% ylabel('PDF','FontSize',12)
% set(gca,'FontSize',12)
% grid on

% Filtered range error pdf - for patent, adapted for a particular data set
windowsize = 10;
n1 = filter(ones(1,windowsize)/windowsize,1,N);
figure(30)
clf
plot(errrange,n1)
hold on
xmin = -2; xmax = 3; ymin = 0.0001; ymax = 0.06; npts = 100;
plot(xmin*ones(1,npts),linspace(ymin,ymax,npts),'r')
plot(linspace(xmin,xmax,npts),ymax*ones(1,npts),'r')
plot(xmax*ones(1,npts),linspace(ymin,ymax,npts),'r')
plot(linspace(minerr,xmin,npts),ymin*ones(1,npts),'r')
plot(linspace(xmax,maxerr,npts),ymin*ones(1,npts),'r')
plot(linspace(-2,0,100), linspace(0.0001,0.06,100),'k')
plot(linspace(0,3,100), linspace(0.06,0.0001,100),'k')

axis([-6 6 0 0.1])
xlabel('Range error (m)','FontSize',12)
ylabel('PDF','FontSize',12)
set(gca,'FontSize',12)
grid on



% Unfiltered range error pdf - histogram
figure(31)
clf
plot(errrange, N);
%hold on
%plot(zeros(size(0:0.001:(max(N/sum(N))+0.005))), 0:0.001:(max(N/sum(N))+0.005), ':r')
xlabel('Range error (m)','FontSize',12)
ylabel('PDF','FontSize',12)
grid on

% CDF of range error distribution
rerr1 = abs(rerr1);
maxerr = max(rerr1);
minerr = 0;
errrange = minerr:(maxerr-minerr)/500:maxerr;
N = histc(rerr1, errrange);
figure(35)
clf
plot(errrange,cumsum(N)/sum(N),'r','LineWidth',2)
xlabel('Range error (m)','FontSize',12)
ylabel('CDF','FontSize',12)
grid on

% probability mass of range error < 1 m  and range error > 1
fprintf('\n Median error is %12.2f\n\n', median(abs(rerr1)));
fprintf('\n Range error <  1 is %d/%d as a percentage this is %12.2f\n', sum(abs(rerr1)<1),length(rerr1),sum(abs(rerr1)<1)/length(rerr1));
fprintf('\n Range error <  5 is %d/%d as a percentage this is %12.2f\n', sum(abs(rerr1)<5),length(rerr1),sum(abs(rerr1)<5)/length(rerr1));
fprintf('\n Range error < 20 is %d/%d as a percentage this is %12.2f\n', sum(abs(rerr1)<20),length(rerr1),sum(abs(rerr1)<10)/length(rerr1));

rerr_ret = rerr(:);


% plot range error distribution between a pair of nodes
node1 = 2;
node2 = 5;
RangeBetweenPairs = [squeeze(RangeData(node1,node2,:)) ; squeeze(RangeData(node2,node1,:))];
RangeBetweenPairs(RangeBetweenPairs==0) = [];
TrueRangeBetweenPairs = dist1(node1,node2);
RangeErrBetweenPairs = RangeBetweenPairs - TrueRangeBetweenPairs;
minerr = min(RangeErrBetweenPairs);
maxerr = max(RangeErrBetweenPairs);
errrange = minerr:(maxerr-minerr)/100:maxerr;
N = histc(RangeErrBetweenPairs, errrange);
figure(350)
clf
plot(errrange,N/sum(N),'r','LineWidth',2)
xlabel('Range error (m)','FontSize',12)
ylabel('CDF','FontSize',12)
grid on

% % find the number of valid ranges in each superframe, will be useful to
% % decide the average connectivity
% nValidSF = size(rerr,2);
% numValidRanges = zeros(1, nValidSF);
% for SFidx = 1:nValidSF
%     
%     trerr = rerr(:,SFidx);
%     trerr(trerr==Inf) = [];
%     
%     numValidRanges(1,SFidx) = length(trerr);
% end
% figure(30)
% clf
% plot(numValidRanges)
