%% Development Information
% MAE 466 Spacecraft Dynamics 
% 
% Function Euler_to_DCM1_Omer is the single axis rotation in the 1st axis
% direction
% 
% input[theta]: theta is the euler angle in radians
%
% output [R1]: is the single rotation axis about the 1st axis
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
 function [R1] = Euler_to_DCM1_Omer(theta)
 R1 = [1, 0, 0; 
       0, cos(theta), sin(theta); 
       0, -sin(theta), cos(theta)];
 end