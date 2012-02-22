function OutVAR = LocalVar(nVAR, InData)
[nS nD] = size(InData);

% moving average of acceleration to get velocity

OutVAR = zeros(nS,nD);

nMAe = ceil(nVAR/2);
nMAb = nMAe - nVAR;

%% Startup code of MA kernel
% Preloads the data from the end of the kernel

XSum = squeeze(zeros(1,nD));
X2Sum = squeeze(zeros(1,nD));
nSum = 0;
if nMAe > 1
    for iS = 1:nMAe-1       % Preload MA filter
        if isnan(InData(iS,1))
        else
          XSum = XSum + InData(iS,:);
          X2Sum = X2Sum + InData(iS,:).^2;
          nSum = nSum + 1;
        end
    end % for
end % if


%% Main Variance code and shutdown code
% If available it adds data at end of kernal, 
% then subtracts obsolete data from before the kernel

for iS = 1:nS
    iB = iS + nMAb-1;
    iE = iS + nMAe-1;
    if iE <= nS
        if isnan(InData(iE,1))
        else
            XSum = XSum + InData(iE,:);
            X2Sum = X2Sum + InData(iE,:).^2;
          nSum = nSum + 1;
        end
    end % if
    if iB > 0
        if isnan(InData(iB,1))
        else
          XSum = XSum - InData(iB,:);
          X2Sum = X2Sum - InData(iB,:).^2;
          nSum = nSum - 1;
        end
    end % if
    OutVAR(iS,:) = X2Sum/nSum - (XSum/nSum).^2;
end % for

end
