% Standard iterative least squares algorithm 
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


function [X, resid] = ILS(Parameters, NumGoodBS, GoodBS, GoodRange, X1)

resid = 1;

Dim = Parameters.SpaceDimension;

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

resid = sqrt(sum(ErrVec.^2));


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

