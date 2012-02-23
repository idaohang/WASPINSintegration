
% Calculates and plots the relative separation error of two independentely tracked WASP
% / INS trials that have a known fixed separation. The wsd files to use are
% specified near the start of the script.

clear all;
offset = 786.3116; % worked out in calcTxTimes.m
true_separation = 1.145;

%main( 'C:\Users\AND48V\Documents\MATLAB\WASP_INS_data\20111220\', 'w220-20111220-rotateinplace-00003.wsd' );
%main( 'C:\Users\AND48V\Documents\MATLAB\WASP_INS_data\20111220\', 'w220-20111220-circle-00005.wsd' );
main( 'C:\Users\AND48V\Documents\MATLAB\WASP_INS_data\20111220\', 'w220-20111220-Lpath-00006.wsd' );

global r Cn2b timestamps_IMU timestamps_WASP WASPloc
r1 = r;
Cn2b1 = Cn2b;
timestamps1 = timestamps_IMU;
WASPloc1 = WASPloc;
timestamps_WASP1 = timestamps_WASP;

%main( 'C:\Users\AND48V\Documents\MATLAB\WASP_INS_data\20111220\', 'w315-20111220-rotateinplace-00003.wsd' );
%main( 'C:\Users\AND48V\Documents\MATLAB\WASP_INS_data\20111220\', 'w315-20111220-circle-00005.wsd' );
main( 'C:\Users\AND48V\Documents\MATLAB\WASP_INS_data\20111220\', 'w315-20111220-Lpath-00006.wsd' );

r2 = r;
Cn2b2 = Cn2b;
timestamps2 = timestamps_IMU;
WASPloc2 = WASPloc;
timestamps_WASP2 = timestamps_WASP;

timestamps1 = timestamps1 - offset;
timestamps_WASP1 = timestamps_WASP1 - offset;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synchronise the WASP positions

% Synchronise the start
WASP1delete = 1:find(timestamps_WASP1(1,:) < timestamps_WASP2(1,1),1,'last');
WASP2delete = 1:find(timestamps_WASP2(1,:) < timestamps_WASP1(1,1),1,'last');

WASPloc1(: , WASP1delete) = [];
WASPloc2(: , WASP2delete) = [];

% Synchronise the end
WASP1delete = find(timestamps_WASP1(1,:) > timestamps_WASP2(1,end),1,'first');
WASP2delete = find(timestamps_WASP2(1,:) > timestamps_WASP1(1,end),1,'first');

WASPloc1(: , WASP1delete:end) = [];
WASPloc2(: , WASP2delete:end) = [];

clear WASP1delete;
clear WASP2delete;

% Now align the data since one or both might be missing observations
if size(timestamps_WASP1,2) < size(timestamps_WASP2,2)
    for t=1:size(timestamps_WASP1,2)
        if( t+1 <= size(timestamps_WASP2,2) &&... 
                abs(timestamps_WASP1(1,t) - timestamps_WASP2(1,t+1)) < abs(timestamps_WASP1(1,t) - timestamps_WASP2(1,t)))
            timestamps_WASP2(:,t)=[];
            WASPloc2(:,t)=[];
        end;
    end;
elseif size(timestamps_WASP2,2) < size(timestamps_WASP1,2)
    for t=1:size(timestamps_WASP2,2)
        if( t+1 <= size(timestamps_WASP1,2) &&...
                abs(timestamps_WASP2(1,t) - timestamps_WASP1(1,t+1)) < abs(timestamps_WASP2(1,t) - timestamps_WASP1(1,t)))
            timestamps_WASP1(:,t)=[];
            WASPloc1(:,t)=[];
        end;
    end;
end;

timestamps_WASP2(:,size(WASPloc1,2)+1:end)=[];
WASPloc2(:,size(WASPloc1,2)+1:end)=[];

timestamps_WASP1(:,size(WASPloc2,2)+1:end)=[];
WASPloc1(:,size(WASPloc2,2)+1:end)=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synchronise the IMU positions

% Synchronise the start
r1delete = 1:find(timestamps1(1,:) < timestamps2(1,1),1,'last');
r2delete = 1:find(timestamps2(1,:) < timestamps1(1,1),1,'last');

r1(: , r1delete) = [];
Cn2b1(: , :, r1delete) = [];
timestamps1(: , r1delete) = [];

r2(: , r2delete) = [];
Cn2b2(: , :, r2delete) = [];
timestamps2(: , r2delete) = [];

% Synchronise the end
r1delete = find(timestamps1(1,:) > timestamps2(1,end),1,'first');
r2delete = find(timestamps2(1,:) > timestamps1(1,end),1,'first');

r1(: , r1delete:end) = [];
Cn2b1(: , :, r1delete:end) = [];
timestamps1(: , r1delete:end) = [];

r2(: , r2delete:end) = [];
Cn2b2(: , :, r2delete:end) = [];
timestamps2(: , r2delete:end) = [];

clear r1delete;
clear r2delete;

% Now align the data since both imu's have a slightly different update rate
if size(timestamps1,2) < size(timestamps2,2)
    for t=1:size(timestamps1,2)
        if( t+1 <= size(timestamps2,2) &&... 
                abs(timestamps1(1,t) - timestamps2(1,t+1)) < abs(timestamps1(1,t) - timestamps2(1,t)))
            timestamps2(:,t)=[];
            r2(:,t)=[];
            Cn2b2(:,:,t)=[];
        end;
    end;
elseif size(timestamps2,2) < size(timestamps1,2)
    for t=1:size(timestamps2,2)
        if( t+1 <= size(timestamps1,2) &&...
                abs(timestamps2(1,t) - timestamps1(1,t+1)) < abs(timestamps2(1,t) - timestamps1(1,t)))
            timestamps1(:,t)=[];
            r1(:,t)=[];
            Cn2b1(:,:,t)=[];
        end;
    end;
end;

timestamps2(:,size(timestamps1,2)+1:end)=[];
r2(:,size(r1,2)+1:end)=[];
Cn2b2(:,:,size(Cn2b1,3)+1:end)=[];

timestamps1(:,size(timestamps2,2)+1:end)=[];
r1(:,size(r2,2)+1:end)=[];
Cn2b1(:,:,size(Cn2b2,3)+1:end)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find error in position and orientation and plot

% Calc separation error for IMU
sep_error = abs(sqrt(sum(abs(r1(1:2,:)-r2(1:2,:)).^2,1)) - true_separation);
mn = mean(sep_error,2);
med = median(sep_error,2);
ccdf = [sort(sep_error); 1:-1/(size(sep_error,2)-1):0]; % complementary cumulative distribution function

% Calc separation error for WASP (no interpolation for timing differences)
sep_error_WASP = abs(sqrt(sum(abs(WASPloc1-WASPloc2).^2,1)) - true_separation);
mn_WASP = mean(sep_error_WASP,2);
med_WASP = median(sep_error_WASP,2);
ccdf_WASP = [sort(sep_error_WASP); 1:-1/(size(sep_error_WASP,2)-1):0]; % complementary cumulative distribution function

% Plot
scrsz = get(0,'ScreenSize'); % Full screen
figure('OuterPosition',[1+scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]);

subplot(3, 1, 1), semilogy(ccdf(1,:),ccdf(2,:)), grid;
hold on;
semilogy(ccdf_WASP(1,:),ccdf_WASP(2,:),'r-'),
xlabel('relative separation error m'),
ylabel('complementary cumulative distribution'),
legend('IMU separation error','WASP separation error');
hold off;

subplot(3, 1, 2), plot(timestamps1 - timestamps1(1), sep_error), grid;
hold on;
plot(timestamps_WASP1 - timestamps_WASP1(1), sep_error_WASP, 'r-'),
text(5,max(sep_error)*0.9,strcat('IMU mean=',num2str(mn),' IMU median=',num2str(med)), 'FontSize', 15),
text(5,max(sep_error)*0.7,strcat('WASP mean=',num2str(mn_WASP),' WASP median=',num2str(med_WASP)), 'FontSize', 15),
xlabel('time in secs'),
ylabel('relative separation error m'),
legend('IMU separation error','WASP separation error');
hold off;

C = zeros(3,3,size(Cn2b1,3));
rpy = zeros(3,size(Cn2b1,3));
for k=1:size(C,3)
    C(:,:,k) = Cn2b1(:,:,k)/Cn2b2(:,:,k);
    rpy(:,k) = rad2deg(DCM2RPY(C(:,:,k)));
end;

subplot(3, 1, 3), plot(timestamps1 - timestamps1(1), rpy), grid
xlabel('time in secs'),
ylabel('orientation difference deg'), legend('R','P','Y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the map

map = loadmap;
set(map, 'OuterPosition', [1 1 scrsz(3)/2 scrsz(4)]);
hold on;

% Reset the axes on the map
xmin = round(min([r1(1,:) r2(1,:)])-1);
xmax = round(max([r1(1,:) r2(1,:)])+1);
ymin = round(min([r1(2,:) r2(2,:)])-1);
ymax = round(max([r1(2,:) r2(2,:)])+1);

if( (ymax-ymin)*1.6 < (xmax-xmin))
    extra = (xmax-xmin)/1.6-(ymax-ymin);
    prev = ylim;
    ymax = min(ymax + extra/2, prev(2));
    ymin = max(ymin - extra/2, prev(1));
else 
    extra = (ymax-ymin)*1.6-(xmax-xmin);
    prev = xlim;
    xmax = min(xmax + extra/2, prev(2));
    xmin = max(xmin - extra/2, prev(1));
end;
axis([xmin; xmax; ymin; ymax]),


jump = 10;
scatter(r1(1,1:jump:end),r1(2,1:jump:end),'g.') % Plot integrated INS position
scatter(r2(1,1:jump:end),r2(2,1:jump:end),'r.') 

% Plot orientation 1
step = 410;
nS = size(Cn2b1,3);
direct = zeros(3,size(1:step:nS,2));
j = 1;
for i=1:step:nS
    direct(:,j) = Cn2b1(:,:,i)'*[1;0;0];
    j = j+1;
end;
quiver(r1(1,1:step:nS),r1(2,1:step:nS),direct(1,:),direct(2,:),'mo','LineWidth',2, 'AutoScaleFactor', 0.1); % show direction

% Plot orientation 2
step = 410;
nS = size(Cn2b2,3);
direct = zeros(3,size(1:step:nS,2));
j = 1;
for i=1:step:nS
    direct(:,j) = Cn2b2(:,:,i)'*[1;0;0];
    j = j+1;
end;
quiver(r2(1,1:step:nS),r2(2,1:step:nS),direct(1,:),direct(2,:),'bo','LineWidth',2, 'AutoScaleFactor', 0.1); % show direction

legend('WASP anchor node', 'node 220 position', 'node 315 position', 'node 220 orientation', 'node 315 orientation');
hold off;

