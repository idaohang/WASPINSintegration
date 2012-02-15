function DCM = RPY2DCM(RPY)
%% Generates a DCM according to given Yaw, Pitch, Roll angles
% Args:
% RPY - Accepts either a column or row vector
% RPY(1) = Roll angle
% RPY(2) = Pitch angle
% RPY(3) = Yaw angle
% Returns:
% DCM - 3x3 Matrix representing rotational transformations to convert 
% a vector in body-frame coordinates into the navigation-frame.
%
% Positive direction of rotation defined as clockwise when viewed
% in the positive direction along the axis
% 
% Euler angle sequence used is (1,2,3), ie. apply the Yaw, Pitch, then
% roll (the aeronautics convenction)

  cY = cos(RPY(3));
  sY = sin(RPY(3));
  cP = cos(RPY(2));
  sP = sin(RPY(2));
  cR = cos(RPY(1));
  sR = sin(RPY(1));
  
  cRcY = cR*cY;
  cRsY = cR*sY;
  sRcY = sR*cY;
  sRsY = sR*sY;
  
% Compute exact DCM rather than product of three matrices to avoid 
% substantial accumulation of rounding errors
  DCM = [(cP*cY)        (cP*sY)           (-sP); ... 
         (sRcY*sP-cRsY) (sRsY*sP+cRcY)  (sR*cP); ... 
         (cRcY*sP+sRsY) (cRsY*sP-sRcY)  (cR*cP)];
     
  DCM = DCM'; % Transpose so its rotating from the b-frame to the n-frame

end
