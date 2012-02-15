function RPY = DCM2RPY(DCM)
%% DCM2YPR: Directional Cosine Matrix to Euler Angle Conversion
% Arguments: 
% DCM - 3x3 Matrix representing rotational matrix that transforms coords
% in the body frame into coords in the navigation frame.
% Returns: 
% RPY - 1x3 (row) Vector containing Yaw, Pitch and Roll angles.
%     YPR(1) - Roll (in radians)
%     YPR(2) - Pitch (in radians)
%     YPR(3) - Yaw (in radians)
% Details:
% Extracts the Euler Angles from a given Direction Cosine Matrix
% Note there are ambiguities when the pitch gets close to 90 deg.
% While these have not been compensated for, this function has only 
% been used for plotting purposes (in order to describe a body's
% orientation in a visually understable way).
% Any attempt to use this function for orientation computation must
% compensate for any issues arising due to singularities occurring
% where pitch is 90 deg.


  DCM = DCM';
  RPY(1) = atan2(DCM(2,3),DCM(3,3));
  RPY(2) = asin(-DCM(1,3));
  RPY(3) = atan2(DCM(1,2),DCM(1,1));
end