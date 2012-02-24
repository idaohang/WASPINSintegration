% Calculate the CRLB for the TOA measurements.
%
% Input:
%   X      - Location of the mobile node
%   BSLocn - Locations of the anchor nodes
%   dim    - Dimension of the space (2 or 3)

% Output:
%   CRLB   - Position error bound.

function CRLB = Calc_CRLB(X, BSLocn, dim)

NumBS = size(BSLocn,1);

DiffLoc = (repmat(X', NumBS,1) - BSLocn)';
DistBS = sqrt(sum(DiffLoc.^2));
H = DiffLoc ./ repmat(DistBS,dim,1);
CRLB = inv(H*H');