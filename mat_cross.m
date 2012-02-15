function [ WQ ] = mat_cross( WG )
%MAT_CROSS Returns the skew-symmetric matrix form of a single 3D vector
%   Also called the cross-product operator
w1 = WG(1);
w2 = WG(2);
w3 = WG(3);

WQ = [  0 -w3  w2; 
       w3   0 -w1;  
      -w2  w1   0];

end

