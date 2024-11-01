%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% =================================================================== %
%                              Part 1.1                               %
% =================================================================== %
   % Noiseless nonlinear simulation
    
    clear; clc; close all;
    addpath('toolbox');
    
    PLOTFLAG = false; % decide whether to produce plots or not

   % Load data and parameters
    load('orbitdetermination-finalproj_data_2023_11_14.mat');
    Nlmks = size(pos_lmks_A,2);
    simParameters;
    
   % Configure integration time and tolerances
    tspan = t0:dt_int:tf_int;
    tol=1.e-14;
    OPTIONS = odeset('RelTol',3*tol,'AbsTol',tol);
    
   % Numerically integrate from initial condition to get state history
    x0 = [r0; rdot0];
    processNoise = zeros(3,1);
    [t_nom,x_nom] = ode45(@(t,x) satDynamicModel(t,x,[],processNoise),tspan,x0,OPTIONS);
    
    
   % Perfect measurement derivation
    y_nom_tbl = [];
    y_full = []; 
    obsvTimes = t0:dt_obs:tf_obs;
    for t_k = obsvTimes % For each observation timestep
        
       % Recover s/c position vector, define khatC vector (cam nadir pointing)
        satpos_k_N = x_nom((t_nom==t_k)~=0,1:3)';    % satellite position in inertial frame
        
       % Get rotation matrix R_CtoN
        R_CtoN_k = R_CtoN(:,:,(t_k/dt_obs)+1);
        ihatC_k_N = R_CtoN_k(:,1); jhatC_k_N = R_CtoN_k(:,2); khatC_k_N = R_CtoN_k(:,3);
    
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
            if u_i >= 0 && u_i <= umax && v_i >= 0 && v_i <= vmax && (disttolmk_k_N'*khatC_k_N) > 0
               % Check whether landmark is facing satellite
                if (lmkipos_k_N'*khatC_k_N) < 0
                   % If so, append [timestamp  lmk_id  u   v] to y_table_ideal
                    y_nom_tbl = [y_nom_tbl;...
                             t_k i u_i v_i]; %#ok<*AGROW>   <--- suppresses size change warning
                end
            end
        
            y_full = [y_full; t_k i u_i v_i]; % save all the meaurments sans the boolean condition (Part 2)
        end
    
    end
    
    if PLOTFLAG
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
    
        markers = {'o','+','*','<','x','s','d','v','^','>'};
        
       % Plot the landmark u pixel measurements
        u_fig = figure;
        for i = 1:10
            % Find rows where the second column is equal to i
            landmark_indices = find(y_nom_tbl(:, 2) == i);
            landmark_id_string = sprintf('Lmk #%d', i);
        
            time_vis = y_nom_tbl(landmark_indices, 1);
            u_loc = y_nom_tbl(landmark_indices, 3);
            plot(time_vis, u_loc, '.', 'DisplayName', landmark_id_string, 'Marker', markers{i});
        
            hold on;
        end
        legend;
        labels(gca,{'Time [s]','u [px]'},'Horizontal Nominal Pixel Position of 10 Landmarks')
        fixfig(u_fig);
        
       % Plot the landmark v pixel measurements
        v_fig = figure; 
        for i = 1:10
            % Find rows where the second column is equal to i
            landmark_indices = find(y_nom_tbl(:, 2) == i);
            landmark_id_string = sprintf('Lmk #%d', i);
        
            time_vis = y_nom_tbl(landmark_indices, 1);
            v_loc = y_nom_tbl(landmark_indices, 4);
            plot(time_vis, v_loc, '.', 'DisplayName', landmark_id_string, 'Marker', markers{i});
        
            hold on;
        end
        legend;
        labels(gca,{'Time [s]','v [px]'},'Vertical Nominal Pixel Position of 10 Landmarks')
        fixfig(v_fig);

    end



%% =================================================================== %
%                              Part 1.2                               %
% =================================================================== %

    % See 'linearizedAmat.m' and 'linearizedCmat.m'


%% =================================================================== %
%                              Part 1.3                               %
% =================================================================== %

    % See 'linearizedFmat.m' and 'linearizedHmat.m'


%% =================================================================== %
%                              Part 1.4                               %
% =================================================================== %

   % Simulation of linearized DT dynamics, using initial state perturbation.
    r_pert = 1.e-5 * [1; 1; 1]; % km
    rdot_pert = 1.e-7 * [1; 1; 1]; % km/s
    x0_pert = [r_pert; rdot_pert];
    
   % Set up container to store DT state history
    x_bar = zeros(length(tspan),6);
    x_bar(1,:) = x0_pert';
    
   % Propagate error state
    for i=1:length(tspan)-1    
        Abar_k = linearizedAmat(mu_A, x_nom(i,:)'); % CT
        Fbar_k = linearizedFmat(dt_int, Abar_k);    % ... to DT
        x_bar(i+1,:) = (Fbar_k*x_bar(i,:)')';  % Propagate, save
    end
    
   % Run another ode45 sim, but this time using the perturbed initial condition 
    processNoise = zeros(3,1);
    [t_pert,x_pert] = ode45(@(t,x) satDynamicModel(t,x,[],processNoise),tspan,x0+x0_pert,OPTIONS);
    
    
   % Error plots
    if PLOTFLAG
       % Position error state history
        eposStateHist = figure;
        subplot(3,1,1);
        plot(tspan,x_bar(:,1),'Color',mlc(1),'DisplayName','$x$');
        labels(gca,{'Time [s]','x [km]'},'Position Error State History (DT, Linearized)');
        subplot(3,1,2);
        plot(tspan,x_bar(:,2),'Color',mlc(2),'DisplayName','$y$');
        labels(gca,{'Time [s]','y [km]'},'');
        subplot(3,1,3);
        plot(tspan,x_bar(:,3),'Color',mlc(3),'DisplayName','$z$');
        labels(gca,{'Time [s]','z [km]'},'');
        fixfig(eposStateHist);
        
       % Velocity error state history
        evelStateHist = figure;
        subplot(3,1,1);
        plot(tspan,x_bar(:,4),'Color',mlc(1),'DisplayName','$\dot{x}$');
        labels(gca,{'Time [s]','$\mathrm{\dot{x}}$ [km/s]'},'Velocity Error State History (DT, Linearized)');
        subplot(3,1,2);
        plot(tspan,x_bar(:,5),'Color',mlc(2),'DisplayName','$\dot{y}$');
        labels(gca,{'Time [s]','$\mathrm{\dot{y}}$ [km/s]'},'');
        subplot(3,1,3);
        plot(tspan,x_bar(:,6),'Color',mlc(3),'DisplayName','$\dot{z}$');
        labels(gca,{'Time [s]','$\mathrm{\dot{z}}$ [km/s]'},'');
        fixfig(evelStateHist);


       % Position: [NOM + ERR] & [NL]
        neposStateHist = figure;
        
        subplot(3,1,1);
        plot(tspan,x_bar(:,1)+x_nom(:,1),'Color',mlc(1),'DisplayName','$\delta x+x_{nom}$'); hold on;
        plot(t_pert,x_pert(:,1),'--','Color','k','DisplayName','$x_{NL}$');
        labels(gca,{'Time [s]','x [km]'},'Full-State Linear vs. Non-Linear Position History'); legend('interpreter','latex');
        
        subplot(3,1,2);
        plot(tspan,x_bar(:,2)+x_nom(:,2),'Color',mlc(2),'DisplayName','$\delta y+y_{nom}$'); hold on;
        plot(t_pert,x_pert(:,2),'--','Color','k','DisplayName','$y_{NL}$');
        labels(gca,{'Time [s]','y [km]'},''); legend('interpreter','latex');
        
        subplot(3,1,3);
        plot(tspan,x_bar(:,3)+x_nom(:,3),'Color',mlc(3),'DisplayName','$\delta z+z_{nom}$'); hold on;
        plot(t_pert,x_pert(:,3),'--','Color','k','DisplayName','$z_{NL}$')
        labels(gca,{'Time [s]','z [km]'},''); legend('interpreter','latex');
        
        fixfig(neposStateHist);
        
       % Velocity: [NOM + ERR] & [NL]
        nevelStateHist = figure;
        subplot(3,1,1);
        plot(tspan,x_bar(:,4)+x_nom(:,4),'Color',mlc(1),'DisplayName','$\delta \dot{x}+\dot{x}_{nom}$'); hold on;
        plot(t_pert,x_pert(:,4),'--','Color','k','DisplayName','$\dot{x}_{NL}$');
        labels(gca,{'Time [s]','$\mathrm{\dot{x}}$ [km/s]'},'Full-State Linear vs. Non-Linear Velocity History'); legend('interpreter','latex');
        
        subplot(3,1,2);
        plot(tspan,x_bar(:,5)+x_nom(:,5),'Color',mlc(2),'DisplayName','$\delta \dot{y}+\dot{y}_{nom}$'); hold on;
        plot(t_pert,x_pert(:,5),'--','Color','k','DisplayName','$\dot{y}_{NL}$');
        labels(gca,{'Time [s]','$\mathrm{\dot{y}}$ [km/s]'},''); legend('interpreter','latex');
        
        subplot(3,1,3);
        plot(tspan,x_bar(:,6)+x_nom(:,6),'Color',mlc(3),'DisplayName','$\delta \dot{z}+\dot{z}_{nom}$'); hold on;
        plot(t_pert,x_pert(:,6),'--','Color','k','DisplayName','$\dot{z}_{NL}$');
        labels(gca,{'Time [s]','$\mathrm{\dot{z}}$ [km/s]'},''); legend('interpreter','latex');
        
        fixfig(nevelStateHist);
    end
    
    
   % Simulate linearized measurements
    y_bar = [];

    for t_k=obsvTimes
        
       % Get rotation matrix R_CtoN_k
        R_CtoN_k = R_CtoN(:,:,(t_k/dt_obs)+1);
    
       % Get rotation matrix R_AtoN_k, rotate landmark position vectors
        theta_k = w_A*t_k;
        R_AtoN_k = rotZ(theta_k);
        pos_lmks_N = R_AtoN_k*pos_lmks_A;

       % Get nominal measurements at time k
        y_nom_tbl_k = y_nom_tbl((y_nom_tbl(:,1)==t_k)~=0,:);

       % Recover nominal s/c state vector at time k
        x_nom_k = x_nom((t_nom==t_k)~=0,:)';
        x_bar_k = x_bar((t_pert==t_k)~=0,:)';

       % Calculate CT jacobian, translate to DT
        Cbar_k = linearizedCmat(f, R_CtoN_k, pos_lmks_N, y_nom_tbl_k, x_nom_k);
        Hbar_k = linearizedHmat(Cbar_k);

        y_bar_k = Hbar_k*x_bar_k;
        y_bar_k = reshape(y_bar_k,[2,length(y_bar_k)/2])';

        y_bar = [y_bar; y_bar_k];
        
    end
    
    y_bar_tbl = [y_nom_tbl(:,1:2) y_bar];
    

    if PLOTFLAG

        % Plot the landmark u pixel measurements
        u_fig = figure;
        for i = 1
            % Find rows where the second column is equal to i
            landmark_indices = find(y_bar_tbl(:, 2) == i);
            landmark_id_string = sprintf('Lmk #%d', i);
        
            time_vis = y_bar_tbl(landmark_indices, 1);
            u_loc = y_bar_tbl(landmark_indices, 3);
            plot(time_vis, u_loc, '.', 'DisplayName', landmark_id_string, 'Marker', markers{i});
        
            hold on;
        end
        legend;
        labels(gca,{'Time [s]','u [px]'},'Horizontal Pixel Position Error of Landmark \#1')
        fixfig(u_fig);
        set(gca,'XTick',0:10*3600:80*3600,'XLim',[0 80*3600]);
    
       % Plot the landmark v pixel measurements
        v_fig = figure; 
        for i = 1
            % Find rows where the second column is equal to i
            landmark_indices = find(y_bar_tbl(:, 2) == i);
            landmark_id_string = sprintf('Lmk #%d', i);
        
            time_vis = y_bar_tbl(landmark_indices, 1);
            v_loc = y_bar_tbl(landmark_indices, 4);
            plot(time_vis, v_loc, '.', 'DisplayName', landmark_id_string, 'Marker', markers{i});
        
            hold on;
        end
        legend;
        labels(gca,{'Time [s]','v [px]'},'Vertical Pixel Position Error of Landmark \#1')
        fixfig(v_fig);
        set(gca,'XTick',0:10*3600:80*3600,'XLim',[0 80*3600],'YTick',-200:50:0,'YLim',[-200 0]);
    
    end
    
    
    
    

