%% Description 
%
% ASEN 5044 Asteroid Project Part 2 
% Use the NEES statistical test to If the KF works properly as 
% per our DT state space model and noise specs 
% (i.e. if it meets the consistency criteria #1-#3 laid out earlier), then we must have:
%
%% Ouput 
%  
% Normalized Estimation Error Sqaured (NEES) Statistic
%
%% Input 
%
% state error, covariance of innovation, and filter covariance
%
%% Version 
% creation date: 4 December 2023
% Matlab version R2023b trial 

%% Revision 
% 
% V1.0 | 4 Dec 2023 | Carlos F Chavez | creation 
% 
%% Program 

function NEESandNIS = calculateNEESNIS(stateError, measurementError, filterCovariance, innovation, alpha, N_MC_Runs, numStates, numMeas)
    % Calculate NEES statistic scalar
    NEES = stateError' / filterCovariance * stateError;

    % Calculate the NIS statistic scalar 
    NIS = measurementError' / innovation * measurementError;

    % Calculate bounds for the NEES 
    r1_NEES = chi2inv(alpha, N_MC_Runs*numStates) / N_MC_Runs;
    r2_NEES = chi2inv(1 - (alpha / 2), N_MC_Runs) / N_MC_Runs;

    % Calculate bounds for the NIS
    r1_NIS = chi2inv(alpha, N_MC_Runs) / N_MC_Runs;
    r2_NIS = chi2inv(1 - (alpha / 2), N_MC_Runs*numMeas) / N_MC_Runs;

    % Plot NEES and bounds
    figure;
    plot(NEES, 'b.-', 'DisplayName', 'NEES');
    hold on;
    plot(r1_NEES, 'r--', 'DisplayName', 'Lower Bound');
    plot(r2_NEES, 'r--', 'DisplayName', 'Upper Bound');
    title('NEES Estimation Results');
    xlabel('Time Step');
    ylabel('NEES');
    legend('show');
    grid on;

    % Plot NIS and bounds
    figure;
    plot(NEES, 'b.-', 'DisplayName', 'NEES');
    hold on;
    plot(r1_NIS, 'r--', 'DisplayName', 'Lower Bound');
    plot(r2_NIS, 'r--', 'DisplayName', 'Upper Bound');
    title('NIS Estimation Results');
    xlabel('Time Step');
    ylabel('NIS');
    legend('show');
    grid on;
end
