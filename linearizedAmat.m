% =================================================================================================
% Function: linearizedAmat
% Calculates the CT state Jacobian at time k.
% =================================================================================================
% Inputs:
%         mu_A   gravitational parameter of the asteroid [km3/s2]
%      x_nom_k   (6 x 1) nominal state at time k
% Outputs:
%       Abar_k   (6 x 6) linearized A matrix, corresponding to xbar_dot_k = Abar_k*xbar_k where
%                        xbar_k = x_k - x_nom_k
% =================================================================================================
function Abar_k = linearizedAmat(mu_A, x_nom_k)

    % Only effects from 2BP acceleration, so only need position state
    pos_k = x_nom_k(1:3);
    x1 = pos_k(1); x2 = pos_k(2); x3 = pos_k(3);

    % Define common denominator for bottom left block
    cd = (pos_k'*pos_k)^(3/2); 

    % Create bottom left block of Abar first
    accelTerms = -mu_A/cd*[x2^2 + x3^2     x1*x2           x1*x3;...
                           x1*x2           x1^2 + x3^2     x2*x3;...
                           x1*x3           x2*x3           x1^2+x2^2];

    % Assemble
    Abar_k = [zeros(3)     eye(3);...
              accelTerms   zeros(3)];

end