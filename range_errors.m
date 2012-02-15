


% Find the true range between nodes
% Note that mobile nodes are not excluded here so interpret with care!
% Need to load a matlab WASP file first 

m = size(NodeList,1);
n = size(RangeData,3);

RangeTrue = zeros(m,m);
for i=1:m
    for j = 1:m
        RangeTrue(i,j) = norm(NodeLocn(i,2:4) - NodeLocn(j,2:4));
    end;
end;

RangeError = zeros(m,m,n);
for i=1:n
    RangeError(:,:,i) = RangeTrue - RangeData(:,:,i);
    for mobile = MobileNodeLocn;
        RangeError(mobile,mobile,i) = 0; % Exclude mobile nodes
    end;
end;

% Make a plot with RSS and SNR, and delete outliers
ErrorPlot = [reshape(RangeError,1,m*m*n); reshape(RSSvalues,1,m*m*n);reshape(rxSNR,1,m*m*n)];
delete = [];
for i=1:m*m*n
    if (abs(ErrorPlot(1,i)) > 1 || ErrorPlot(1,i)==0)
        delete = [delete i];
    end;
end;
ErrorPlot(:,delete) = [];
clear delete;

% Some statistics
meanRangeErr = mean(ErrorPlot(1,:))
varRangeErr = var(ErrorPlot(1,:))


% Make the plots - first a histogram
figure;
x = -1:0.05:1;
hist(ErrorPlot(1,:),x);
title('WASP range error distribution (m)');

figure;
scattercloud(ErrorPlot(1,:),ErrorPlot(2,:));

figure;
plot(ErrorPlot(1,:),ErrorPlot(2,:),'r+'),grid;
ylabel('RSSvalue'),
xlabel('WASP range error (m)');

figure;
plot(ErrorPlot(1,:),ErrorPlot(3,:),'b+'),grid;
ylabel('rxSNR'),
xlabel('WASP range error (m)');

figure;
plot(ErrorPlot(1,:).^2,ErrorPlot(2,:),'r+'),grid;
ylabel('RSSvalue'),
xlabel('WASP range squared error (m)');