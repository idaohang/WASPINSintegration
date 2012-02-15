function OutMA = MovAvg2(nMA, InData)
[nS nD] = size(InData);

% moving average of acceleration to get velocity

OutMA = zeros(nS,nD);

nMAe = ceil(nMA/2);
nMAb = nMAe - nMA;

%% Startup code of MA kernel
% Preloads the data from the end of the kernel

ASum = squeeze(zeros(1,nD));
nSum = 0;
if nMAe > 1
    for iS = 1:nMAe-1       % Preload MA filter
        if isnan(InData(iS,1))
        else
          ASum = ASum + InData(iS,:);
          nSum = nSum + 1;
%          x = 0;
%          display(horzcat(x, iS, nSum, ASum));
        end
    end % for
end % if

%logMA = zeros(3,ss1);

%% Main MA code and shutdown code
% If available it adds data at end of kernal, 
% then subtracts obsolete data from before the kernel

for iS = 1:nS
    iB = iS + nMAb-1;
    iE = iS + nMAe-1;
    if iE <= nS
        if isnan(InData(iE,1))
        else
          ASum = ASum + InData(iE,:);
          nSum = nSum + 1;
%          x = 1;
%          display(horzcat(x, iS, iE, nSum, ASum));
        end
    end % if
    if iB > 0
        if isnan(InData(iB,1))
        else
          ASum = ASum - InData(iB,:);
          nSum = nSum - 1;
%          x = 2;
%          display(horzcat(x, iS, iB, nSum, ASum));
        end
    end % if
    OutMA(iS,:) = ASum/nSum;
end % for

end
