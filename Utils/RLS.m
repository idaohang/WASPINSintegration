% Robust least squares algorithm.
%
% Input:
%   Parameters - A structure consisting of various parameters.
%   NumGoodBS  - Number of anchor nodes.
%   GoodBS     - Locations of the anchor nodes.
%   GoodRange  - Range measurement to each anchor node
%   X1         - Initial guess to start the iteration. If it is empty, the
%                a linear least squares is performed
%
% Output:
%   X          - Estimated location.
%   resid      - Residual of the location estiamte
%
%
function [X, resid] = RLS(Parameters, NumGoodBS, GoodBS, GoodRange, X1)

resid = -1;

Dim = Parameters.SpaceDimension;
GoodPosnError = Parameters.RLS_GoodPosnError;  % if the estimated position error is less than this don't seek better solution
MinBS = Parameters.MinimumNumberOfAnchors;	   % minimum number of base stations for solution
RemoveMargin = Parameters.RLS_RemoveMargin;	   % to remove a base station the estimated postion error must be improved by this factor

if (NumGoodBS < Parameters.MinimumNumberOfAnchors)
    X = [];
    return;
end

TrialGoodRange = GoodRange;
TrialGoodBS = GoodBS;
TrialNumGoodBS = NumGoodBS;

% calculate linear solution
if nargin==4
    X = LinearSolution(TrialNumGoodBS, TrialGoodBS, TrialGoodRange, Dim);
elseif nargin==5
    X = X1;
end

% non-linear refinement
[X,ErrVec] = NonLinearSolution(X, TrialNumGoodBS, TrialGoodBS, TrialGoodRange, Dim);
% estimate location error
PosnError = CRLB_TOA(X, TrialNumGoodBS, TrialGoodBS, Dim)*sqrt(sum(ErrVec.^2)/(length(ErrVec)-2));

% iteratively remove basestations until no further improvement
while (PosnError > GoodPosnError) && (TrialNumGoodBS > MinBS)

    BSbest = 0;

    % take each BS in turn, and keep solution with lowest error
    for BSremove=1:NumGoodBS
        
        % temporary variables with base station removed
        GBSmask = ones(1,NumGoodBS);
        GBSmask(BSremove) = 0;
        GBSidx = find(GBSmask==1);
        TrialGoodRange = GoodRange(GBSidx);
        TrialGoodBS = GoodBS(GBSidx,:);
        TrialNumGoodBS = NumGoodBS - 1;
        
        % find solution and estimated position error
        [X2,ErrVec2] = NonLinearSolution(X, TrialNumGoodBS, TrialGoodBS, TrialGoodRange, Dim);
        TrialPosnError = CRLB_TOA(X2, TrialNumGoodBS, TrialGoodBS, Dim)*sqrt(sum(ErrVec2.^2)/(length(ErrVec2)-2));
        
        % keep best result
        if TrialPosnError*RemoveMargin < PosnError
            BSbest = BSremove;
            PosnError = TrialPosnError;
            X = X2;
            BestGoodRange = TrialGoodRange;
            BestGoodBS = TrialGoodBS;
            BestNumGoodBS = TrialNumGoodBS;
            ErrVec = ErrVec2;
        end
    end
    
    % exit from loop if no improvement
    if BSbest == 0
        break;
    else
        GoodRange = BestGoodRange;
        GoodBS = BestGoodBS;
        NumGoodBS = BestNumGoodBS;
    end
end

resid = sum(ErrVec.^2)/length(ErrVec);

% calculate Cramer Rao lower bound (GDOP) for TOA measurement (i.e. fixed delay)
function GDOP = CRLB_TOA(X, TrialNumGoodBS, TrialGoodBS, Dim)
    
%DiffLoc = (repmat(X(1:2)',TrialNumGoodBS,1) - TrialGoodBS)';
DiffLoc = (X(:,ones(TrialNumGoodBS,1))' - TrialGoodBS)';
DistBS = sqrt(sum(DiffLoc.^2));
H = DiffLoc ./ DistBS(ones(Dim,1),:);
GDOP = sqrt(trace(inv(H*H')));

% calculate linear approximation for location
function X = LinearSolution(TrialNumGoodBS, TrialGoodBS, TrialGoodRange, Dim)

AP = 2*( TrialGoodBS(2:end,:) - repmat(TrialGoodBS(1,:),TrialNumGoodBS-1,1) );
g = -TrialGoodRange(2:end).^2 + TrialGoodRange(1)^2 + sum(TrialGoodBS(2:end,:).^2,2) - sum(TrialGoodBS(1,:).^2);
AIV = AP'*AP;

DetAIV = det(AIV);
if DetAIV > 1e-6
    InvAIV = inv(AIV);
    X = InvAIV * AP' * g;
else
    X = zeros(Dim,1);
end
 

% non-linear iterations to refine location. Uses Taylor expansion.
function [X,ErrVec] = NonLinearSolution(X, TrialNumGoodBS, TrialGoodBS, TrialGoodRange, Dim)

for k = 1:5			% five iterations is more than sufficient
    Pdiff = TrialGoodBS - X(:,ones(TrialNumGoodBS,1))';
    Pdist = sqrt(sum(Pdiff.^2,2));
    Mdiff = TrialGoodRange - Pdist;
    H = Pdiff ./ repmat(Pdist,1,Dim);
    AIV = H'*H;
    %DetAIV = AIV(1,1)*AIV(2,2) - AIV(1,2)*AIV(2,1);
    DetAIV = det(AIV);
    if DetAIV > 1e-8
        %InvAIV = [AIV(2,2) -AIV(1,2); -AIV(2,1) AIV(1,1)]/DetAIV;
        InvAIV = inv(AIV);
        dX = InvAIV*H'*Mdiff;
    else
        disp('Warning: Taylor iteration ill conditioned');
        break;
    end
    if sum(dX.^2) < 1e-6
        break;
    end
    if sum(dX.^2) > 1e4
        disp('Warning: Non-linear refinement diverging');
        break;
    end
    X = X - dX;
end

Pdiff = TrialGoodBS - X(:,ones(TrialNumGoodBS,1))';
Pdist = sqrt(sum(Pdiff.^2,2));
ErrVec = TrialGoodRange - Pdist;

