% Resampling algorithm(s) used in the particle filter.
%
% Input:
%   w       - Normalized weights of the particles.
%   Alg     - Chosen resampling algorithm
%
% Output:
%   keepids - IDs of the particles to keep.
%
%

function keepids = ResamplingAlgorithms(w,Alg)

if strcmp(Alg,'SystematicResampling')
    Ns = length(w);

    weight_cdf = [0 cumsum(w)];
    %weight_cdf(end) = [];
    
    k1 = 1;
    u1 = rand/Ns;
    u = u1;
    
    keepids = zeros(Ns,1);
    
    for k2 = 1:Ns
        while u > weight_cdf(k1)
            k1 = k1 + 1;
        end
        keepids(k2) = k1 - 1; %k1 in the tutorial
        u = u1 + (k2-1)/Ns;
    end
else
    error('Unknown resampling method.');
end
