%% Development Information
% MAE 466 Spacecraft Dynamics 
% 
% Function Euler_to_DCM2_Omer is the single axis rotation in the 2nd axis
% direction
% 
% input[theta]: theta is the euler angle in radians
%
% output [R2]: is the single rotation axis about the 2nd axis
% 
% 
% Primary Developer Contact Information:
% Mazin Omer
% Mechanical and Aerospace Engineering '26
% mo00024@mix.wvu.edu
% 5073198922
%
%
% Development History
% Date            Developer        Comments
% -------------   -------------    --------------------------------
% Sep.23rd 2025       M. Omer
 function [R2] = Euler_to_DCM2_Omer(theta)
 R2 = [cos(theta), 0, -sin(theta);
          0,        1, 0;
         sin(theta), 0, cos(theta)];
 end