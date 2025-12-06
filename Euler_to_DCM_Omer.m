%% Development Information
% MAE 466 Spacecraft Dynamics 
% 
% Function Euler_to_DCM_Omer finds DCM matrix from Euler angles only for
%                            313 and 321 sequences
% 
% input [theta_vec]:    theta vector is the euler vector in radians
%       [sequence]:     is the number representing the euler angle sequence  
%
% output [R]:           3x3 DCM rotation Matrix 
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
function [R] = Euler_to_DCM_Omer(theta_vec,sequence)

% Initialize the rotation matrix in case we use in simulink
    R = eye(3);

     if strcmpi(sequence, '313') == 1
        R_first  = Euler_to_DCM3_Omer(theta_vec(1));
        R_middle = Euler_to_DCM1_Omer(theta_vec(2));
        R_end    = Euler_to_DCM3_Omer(theta_vec(3));

    elseif strcmpi(sequence, '321') == 1
        R_first  = Euler_to_DCM3_Omer(theta_vec(1));
        R_middle = Euler_to_DCM2_Omer(theta_vec(2));
        R_end    = Euler_to_DCM1_Omer(theta_vec(3));

    else
        error('This function only supports sequence 313 and 321');
    end

    R = R_end * R_middle * R_first;

% replacing small values with 0
   for i= 1:9
       if abs(R(i)) < 1e-10
           R(i) = 0;
       end
   end
    
end