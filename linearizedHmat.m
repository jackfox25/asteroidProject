% =================================================================================================
% Function: linearizedHmat
% Calculates the DT measurement Jacobian at time k.
% =================================================================================================
% Inputs:
%       Cbar_k   (2*Nlmks_k x 6) linearized C matrix, corresponding to ybar_k = Cbar_k*xbar_k where
%                                ybar_k = y_k - y_nom_k , xbar_k = x_k - x_nom_k
% Outputs:
%       Hbar_k   (2*Nlmks_k x 6) linearized H matrix, corresponding to ybar_k = Hbar_k*xbar_k where
%                                ybar_k = y_k - y_nom_k , xbar_k = x_k - x_nom_k
% =================================================================================================
function Hbar_k = linearizedHmat(Cbar_k)

    Hbar_k = Cbar_k;

end