% Returns the motion models. Only nearly constant velocity model is used by
% various processing algorithms.
% 
% Input:
%   First argument  - space dimension
%   Other arguments - A pair consisting of a string (name of the variable)
%                     and sampling time.
%
% Output:
%   Model(s) requested
%

function varargout = GetMotionModel(varargin)

DIM = varargin{1};

paramCnt = 1;
argout = 0;

while paramCnt < nargin
    
    paramCnt = paramCnt + 1;

    if strcmp(varargin{paramCnt},'StateTransitionMatrix_CV')
        paramCnt = paramCnt + 1;
        dT = varargin{paramCnt};
        Phi = [1 dT;
               0  1];
        if DIM == 2
            Phi = blkdiag(Phi,Phi);
        elseif DIM == 3
            Phi = blkdiag(Phi,Phi,Phi);
        end
        
        argout = argout + 1;
        varargout{argout} = Phi;
    elseif strcmp(varargin{paramCnt},'StateTransitionMatrix_CA')
        paramCnt = paramCnt + 1;
        dT = varargin{paramCnt};
        Phi = [1 dT dT^2/2;
               0 1  dT    ;
               0 0  1    ];
        if DIM == 2
            Phi = blkdiag(Phi,Phi);
        elseif DIM == 3
            Phi = blkdiag(Phi,Phi,Phi);
        end
        
        argout = argout + 1;
        varargout{argout} = Phi;
    elseif strcmp(varargin{paramCnt},'ProcessNoiseGain_CV')
        paramCnt = paramCnt + 1;
        dT = varargin{paramCnt};
        TAU = [dT;
                1];
        if DIM == 2
            TAU = blkdiag(TAU,TAU);
        elseif DIM == 3
            TAU = blkdiag(TAU,TAU,TAU);
        end
        
        argout = argout + 1;
        varargout{argout} = TAU;
    elseif strcmp(varargin{paramCnt},'ProcessNoiseGain_CA')
        paramCnt = paramCnt + 1;
        dT = varargin{paramCnt};
        TAU = [(1/2)*dT^2;
                     dT  ;
                      1  ];
        if DIM == 2
            TAU = blkdiag(TAU,TAU);
        elseif DIM == 3
            TAU = blkdiag(TAU,TAU,TAU);
        end
        
        argout = argout + 1;
        varargout{argout} = TAU;
    elseif strcmp(varargin{paramCnt},'ProcessNoiseCovarianceMatrix_CV')
        paramCnt = paramCnt + 1;
        dT = varargin{paramCnt};
        paramCnt = paramCnt + 1;
        q = varargin{paramCnt};        
        TAU = q * [dT^3/3  dT^2/2;
                   dT^2/2  dT];   
        if DIM == 2
            Q = blkdiag(TAU,TAU);
        elseif DIM == 3
            Q = blkdiag(TAU,TAU,TAU);
        end              
        
        argout = argout + 1;
        varargout{argout} = Q;        
    elseif strcmp(varargin{paramCnt},'ProcessNoiseCovarianceMatrix_CA')
        paramCnt = paramCnt + 1;
        dT = varargin{paramCnt};
        paramCnt = paramCnt + 1;
        q = varargin{paramCnt};
        TAU = [(1/2)*dT^2;
                     dT  ;
                      1  ];
        TAU = blkdiag(TAU,TAU);
        Q = q*TAU*TAU';
        
        argout = argout + 1;
        varargout{argout} = Q;        
    end
end

