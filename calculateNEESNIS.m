% =================================================================================================
% Function: calculateNEESNIS
% Calculates the NEES and NIS criteria at time k.
% =================================================================================================
% Inputs:
%      x_bar_k   (6 x 1) state error at time k
%      y_bar_k   (2*Nlmks x 1) measurement error at time k
%         P_kp   (6 x 6) state covariance matrix (a posteriori)
%            S   (xx x xx) innovation covariance matrix ( = H*Pkm*H'+R )
%        alpha   plot threshold (i.e. 0.05 for +/-2sigma)
%          NMC   # of Monte-Carlo runs
% Outputs:
%         NEES   structure containing fields:
%                  - NEES, r1, r2    (i.e. NEES.NEES or NEES.r1)
%          NIS   structure containing fields:
%                  - NIS, r1, r2     (i.e. NIS.NIS or NIS.r1)
% =================================================================================================
function [NEES, NIS] = calculateNEESNIS(x_bar_k, y_bar_k, P_kp, S, alpha, NMC)
    
    n = size(x_bar_k,1); p = size(y_bar_k,1);

    % Calculate NEES statistic scalar
    NEES.NEES = x_bar_k' / P_kp * x_bar_k;

    % Calculate the NIS statistic scalar 
    NIS.NIS = y_bar_k' / S * y_bar_k;

    % Calculate bounds for the NEES 
    NEES.r1 = chi2inv(alpha/2, NMC*n) / NMC;
    NEES.r2 = chi2inv(1 - (alpha/2), NMC*n) / NMC;

    % Calculate bounds for the NIS
    NIS.r1 = chi2inv(alpha/2, NMC*p) / NMC;
    NIS.r2 = chi2inv(1 - (alpha/2), NMC*p) / NMC;
    
end
