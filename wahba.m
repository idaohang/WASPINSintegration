function [ Cb2n ] = wahba( Ab, An )
%wahba  Finds a Direction Cosine Matrix between two coordinate systems, 
%   given a set of paired vector observations (Wahba's problem).
%   Inputs to the function are:
%       two sets of 3 dimensional vectors An and Ab in mxn matrix form,
%       i.e. m = 3, n = number of vectors.
%   Output is the DCM from the b-frame to the n-frame (Ab to An). 

ns = size(An,2);
B = zeros(3,3);
for t=1:ns
   B = B + An(:,t)*Ab(:,t)';
end

[U S V] = svd(B);
M = diag([1; 1; det(U)*det(V')]);
Cb2n = (U*M*V')';

end

