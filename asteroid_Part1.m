%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% =================================================================== %
%                              Part 1.1                               %
% =================================================================== %
    % Noiseless nonlinear simulation
    
    clear; clc; close all;
    addpath('toolbox');
    
    % Load data and parameters
    load('orbitdetermination-finalproj_data_2023_11_14.mat');
    Nlmks = size(pos_lmks_A,2);
    simParameters;
    
    % Configure integration time and tolerances
    tspan = t0:dt_int:tf;
    tol=1.e-14;
    OPTIONS = odeset('RelTol',3*tol,'AbsTol',tol);
    
    % Numerically integrate from initial condition to get state history
    x0 = [r0; rdot0];
    processNoise = zeros(3,1);
    [t_nom,x_nom] = ode45(@(t,x) satDynamicModel(t,x,[],processNoise),tspan,x0,OPTIONS);
    
    
    % Perfect measurement derivation
    y_nom = [];
    
    obsvTimes = unique(y_table(:,1))';
    for t_k = obsvTimes % For each observation timestep
        
        % Recover s/c position vector, define khatC vector (cam nadir pointing)
        satpos_k_N = x_nom(find(t_nom==t_k),1:3)';          % satellite position in inertial frame
        khatC_k_N = -satpos_k_N/norm(satpos_k_N);   % unit vector of camera pointing axis in inertial frame
        
        % Get rotation matrix R_CtoN
        R_CtoN_k = R_CtoN(:,:,(t_k/dt_obs)+1);
        ihatC_k_N = R_CtoN_k(:,1); jhatC_k_N = R_CtoN_k(:,2);
    
        % Get rotation matrix R_AtoN
        theta = w_A*t_k;
        R_AtoN_k = rotZ(theta);
    
        for i = 1:Nlmks % For each landmark
           
           % Get landmark position in intertial frame
            lmkipos_k_A = pos_lmks_A(:,i);
            lmkipos_k_N = R_AtoN_k*lmkipos_k_A;
            disttolmk_k_N = lmkipos_k_N - satpos_k_N;
    
           % Simulate ideal measurement
            u_i = f * (disttolmk_k_N'*ihatC_k_N) / (disttolmk_k_N'*khatC_k_N) + u0;
            v_i = f * (disttolmk_k_N'*jhatC_k_N) / (disttolmk_k_N'*khatC_k_N) + v0;
    
           % Check whether landmark is in camera FOV
            if u_i >= 0 && u_i <= umax && (disttolmk_k_N'*khatC_k_N) > 0
               % Check whether landmark is facing satellite
                if dot(lmkipos_k_N',khatC_k_N) < 0
                   % If so, append [timestamp  lmk_id  u   v] to y_table_ideal
                    y_nom = [y_nom;...
                             t_k i u_i v_i]; %#ok<*AGROW>   <--- suppresses size change warning
                end
            end
    
        end
    
    end
    
    % Position state history
    posStateHist = figure;
    subplot(3,1,1);
    plot(t_nom,x_nom(:,1),'Color',mlc(1),'DisplayName','$x$');
    labels(gca,{'Time [s]','x [km]'},'NL Nominal Position State History');
    subplot(3,1,2);
    plot(t_nom,x_nom(:,2),'Color',mlc(2),'DisplayName','$y$');
    labels(gca,{'Time [s]','y [km]'},'');
    subplot(3,1,3);
    plot(t_nom,x_nom(:,3),'Color',mlc(3),'DisplayName','$z$');
    labels(gca,{'Time [s]','z [km]'},'');
    fixfig(posStateHist);
    
    % Velocity state history
    velStateHist = figure;
    subplot(3,1,1);
    plot(t_nom,x_nom(:,4),'Color',mlc(1),'DisplayName','$\dot{x}$');
    labels(gca,{'Time [s]','$\mathrm{\dot{x}}$ [km/s]'},'NL Nominal Velocity State History');
    subplot(3,1,2);
    plot(t_nom,x_nom(:,5),'Color',mlc(2),'DisplayName','$\dot{y}$');
    labels(gca,{'Time [s]','$\mathrm{\dot{y}}$ [km/s]'},'');
    subplot(3,1,3);
    plot(t_nom,x_nom(:,6),'Color',mlc(3),'DisplayName','$\dot{z}$');
    labels(gca,{'Time [s]','$\mathrm{\dot{z}}$ [km/s]'},'');
    fixfig(velStateHist);
    
    markers = {'o','x','+','*','s','d','v','^','<','>'};
    
    % Plot the landmark u pixel measurements
    u_fig = figure;
    for i = 1:10
        % Find rows where the second column is equal to i
        landmark_indices = find(y_nom(:, 2) == i);
        landmark_id_string = sprintf('Lmk #%d', i);
    
        time_vis = y_nom(landmark_indices, 1);
        u_loc = y_nom(landmark_indices, 3);
        plot(time_vis, u_loc, '.', 'DisplayName', landmark_id_string, 'Marker', markers{i});
    
        hold on;
    end
    legend;
    labels(gca,{'Time [s]','u [px]'},'Horizontal Pixel Position of 10 Landmarks (Ideal Measurements)')
    fixfig(u_fig);
    
    % Plot the landmark v pixel measurements
    v_fig = figure; 
    for i = 1:10
        % Find rows where the second column is equal to i
        landmark_indices = find(y_nom(:, 2) == i);
        landmark_id_string = sprintf('Lmk #%d', i);
    
        time_vis = y_nom(landmark_indices, 1);
        v_loc = y_nom(landmark_indices, 4);
        plot(time_vis, v_loc, '.', 'DisplayName', landmark_id_string, 'Marker', markers{i});
    
        hold on;
    end
    legend;
    labels(gca,{'Time [s]','v [px]'},'Vertical Pixel Position of 10 Landmarks (Ideal Measurements)')
    fixfig(v_fig);





% =================================================================== %
%                              Part 1.2                               %
% =================================================================== %

    % See 'linearizedAmat.m' and 'linearizedCmat.m'


% =================================================================== %
%                              Part 1.3                               %
% =================================================================== %

    % See 'linearizedFmat.m' and 'linearizedGmat.m'





%% =================================================================== %
%                              Part 1.4                               %
% =================================================================== %

    % Simulation of linearized DT dynamics, using initial state perturbation.
    r_pert = 1.e-5 * [1; 1; 1]; % km
    rdot_pert = 1.e-7 * [1; 1; 1]; % km/s
    x0_pert = [r_pert; rdot_pert];
    
    % Set up container to store DT state history
    xbar_hist = zeros(length(tspan),6);
    xbar_hist(1,:) = x0_pert';
    
    for i=1:length(tspan)-1    
        Abar_k = linearizedAmat(mu_A, x_nom(i,:)'); % CT
        Fbar_k = linearizedFmat(dt_int, Abar_k);    % ... to DT
        xbar_hist(i+1,:) = (Fbar_k*xbar_hist(i,:)')';  % Propagate, save
    end
    
    % Run another ode45 sim, but this time using the perturbed initial condition 
    processNoise = zeros(3,1);
    [t_pert,x_pert] = ode45(@(t,x) satDynamicModel(t,x,[],processNoise),tspan,x0+x0_pert,OPTIONS);
    
    
    % Error plots
    if true
        % Position error state history
        eposStateHist = figure;
        subplot(3,1,1);
        plot(tspan,xbar_hist(:,1),'Color',mlc(1),'DisplayName','$x$');
        labels(gca,{'Time [s]','x [km]'},'Position Error State History');
        subplot(3,1,2);
        plot(tspan,xbar_hist(:,2),'Color',mlc(2),'DisplayName','$y$');
        labels(gca,{'Time [s]','y [km]'},'');
        subplot(3,1,3);
        plot(tspan,xbar_hist(:,3),'Color',mlc(3),'DisplayName','$z$');
        labels(gca,{'Time [s]','z [km]'},'');
        fixfig(eposStateHist);
        
        % Velocity error state history
        evelStateHist = figure;
        subplot(3,1,1);
        plot(tspan,xbar_hist(:,4),'Color',mlc(1),'DisplayName','$\dot{x}$');
        labels(gca,{'Time [s]','$\mathrm{\dot{x}}$ [km/s]'},'Velocity Error State History');
        subplot(3,1,2);
        plot(tspan,xbar_hist(:,5),'Color',mlc(2),'DisplayName','$\dot{y}$');
        labels(gca,{'Time [s]','$\mathrm{\dot{y}}$ [km/s]'},'');
        subplot(3,1,3);
        plot(tspan,xbar_hist(:,6),'Color',mlc(3),'DisplayName','$\dot{z}$');
        labels(gca,{'Time [s]','$\mathrm{\dot{z}}$ [km/s]'},'');
        fixfig(evelStateHist);
    end
    
    
    % Position: [NOM + ERR] & [NL]
    neposStateHist = figure;
    
    subplot(3,1,1);
    plot(tspan,xbar_hist(:,1)+x_nom(:,1),'Color',mlc(1),'DisplayName','$\delta x+x_{nom}$'); hold on;
    plot(t_pert,x_pert(:,1),'--','Color','k','DisplayName','$x_{NL}$');
    labels(gca,{'Time [s]','x [km]'},'Position Error State History'); legend('interpreter','latex');
    
    subplot(3,1,2);
    plot(tspan,xbar_hist(:,2)+x_nom(:,2),'Color',mlc(2),'DisplayName','$\delta y+y_{nom}$'); hold on;
    plot(t_pert,x_pert(:,2),'--','Color','k','DisplayName','$y_{NL}$');
    labels(gca,{'Time [s]','y [km]'},''); legend('interpreter','latex');
    
    subplot(3,1,3);
    plot(tspan,xbar_hist(:,3)+x_nom(:,3),'Color',mlc(3),'DisplayName','$\delta z+z_{nom}$'); hold on;
    plot(t_pert,x_pert(:,3),'--','Color','k','DisplayName','$z_{NL}$')
    labels(gca,{'Time [s]','z [km]'},''); legend('interpreter','latex');
    
    fixfig(neposStateHist);
    
    % Velocity: [NOM + ERR] & [NL]
    nevelStateHist = figure;
    subplot(3,1,1);
    plot(tspan,xbar_hist(:,4)+x_nom(:,4),'Color',mlc(1),'DisplayName','$\delta \dot{x}+\dot{x}_{nom}$'); hold on;
    plot(t_pert,x_pert(:,4),'--','Color','k','DisplayName','$\dot{x}_{NL}$');
    labels(gca,{'Time [s]','$\mathrm{\dot{x}}$ [km/s]'},'Velocity Error State History'); legend('interpreter','latex');
    
    subplot(3,1,2);
    plot(tspan,xbar_hist(:,5)+x_nom(:,5),'Color',mlc(2),'DisplayName','$\delta \dot{y}+\dot{y}_{nom}$'); hold on;
    plot(t_pert,x_pert(:,5),'--','Color','k','DisplayName','$\dot{y}_{NL}$');
    labels(gca,{'Time [s]','$\mathrm{\dot{y}}$ [km/s]'},''); legend('interpreter','latex');
    
    subplot(3,1,3);
    plot(tspan,xbar_hist(:,6)+x_nom(:,6),'Color',mlc(3),'DisplayName','$\delta \dot{z}+\dot{z}_{nom}$'); hold on;
    plot(t_pert,x_pert(:,6),'--','Color','k','DisplayName','$\dot{z}_{NL}$');
    labels(gca,{'Time [s]','$\mathrm{\dot{z}}$ [km/s]'},''); legend('interpreter','latex');
    
    fixfig(nevelStateHist);
    
    
    % Simulate linearized measurements
    
    % TODO
    
    
    
    
    
    
    
    
    

