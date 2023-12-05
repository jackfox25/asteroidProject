% =================================================================================================
% Function: linearizedFmat
% Calculates the DT state Jacobian at time k.
% =================================================================================================
% Inputs:
%           dt   state integration timestep [s]         
%       Abar_k   (6 x 6) linearized A matrix, corresponding to xbar_dot_k = Abar_k*xbar_k where
%                        xbar_k = x_k - x_nom_k
% Outputs:
%       Fbar_k   (6 x 6) linearized F matrix, corresponding to xbar_kp1 = Fbar_k*xbar_k
% =================================================================================================
function Fbar_k = linearizedFmat(dt, Abar_k)

    Fbar_k = eye(6) + dt*Abar_k;

end