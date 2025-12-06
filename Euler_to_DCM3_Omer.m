%% Development Information
% MAE 466 Spacecraft Dynamics 
% 
% Function Euler_to_DCM3_Omer is the single axis rotation in the 3rd axis
% direction
% 
% input[theta]: theta is the euler angle in radians
%
% output [R3]: is the single rotation axis about the 3nd axis
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
 function [R3] = Euler_to_DCM3_Omer(theta)
 R3 = [cos(theta), sin(theta), 0; 
      -sin(theta),  cos(theta), 0; 
       0,           0,          1];
end
 