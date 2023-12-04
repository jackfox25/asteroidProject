%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% =================================================================== %
% %                            Part 1.1                             % %
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
[t,X] = ode45(@(t,x) satDynamicModel(t,x,[],processNoise),tspan,x0,OPTIONS);


% Perfect measurement derivation
y_nom = [];

obsvTimes = unique(y_table(:,1));
for t_k = obsvTimes' % For each observation timestep
    
    % Recover s/c position vector, define khatC vector (cam nadir pointing)
    satpos_k_N = X(find(t==t_k),1:3)';          % satellite position in inertial frame
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
plot(t,X(:,1),'Color',mlc(1),'DisplayName','$x$');
labels(gca,{'Time [s]','x [km]'},'Position State History');
subplot(3,1,2);
plot(t,X(:,2),'Color',mlc(2),'DisplayName','$y$');
labels(gca,{'Time [s]','y [km]'},'');
subplot(3,1,3);
plot(t,X(:,3),'Color',mlc(3),'DisplayName','$z$');
labels(gca,{'Time [s]','z [km]'},'');
fixfig(posStateHist);

% Velocity state history
velStateHist = figure;
subplot(3,1,1);
plot(t,X(:,4),'Color',mlc(1),'DisplayName','$\dot{x}$');
labels(gca,{'Time [s]','$\mathrm{\dot{x}}$ [km/s]'},'Velocity State History');
subplot(3,1,2);
plot(t,X(:,5),'Color',mlc(2),'DisplayName','$\dot{y}$');
labels(gca,{'Time [s]','$\mathrm{\dot{y}}$ [km/s]'},'');
subplot(3,1,3);
plot(t,X(:,6),'Color',mlc(3),'DisplayName','$\dot{z}$');
labels(gca,{'Time [s]','$\mathrm{\dot{z}}$ [km/s]'},'');
fixfig(velStateHist);

markers = {'o','x','+','*','s','d','v','^','<','>'};
% Plot the Landmark u Pixel Measurements
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

% Plot the Landmark v Pixel Measurements
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
% %                            Part 1.2                             % %
% =================================================================== %








