%% PROBLEM 1 %%
% Noiseless nonlinear simulation

clear; clc; close all;
ANIMATE_FLAG = false;

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
y_table_ideal = [];

obsvTimes = unique(y_table(:,1));
for t_k = obsvTimes' % For each observation timestep
    
    % Recover s/c position vector, define khatC vector (cam nadir pointing)
    satpos_k_N = X(find(t==t_k),1:3)';          % satellite position in inertial frame
    khatC_k_N = -satpos_k_N/norm(satpos_k_N);   % unit vector of camera pointing axis in inertial frame
    
    % Get rotation matrix R_CtoN
    R_CtoN_k = R_CtoN(:,:,(t_k/dt_obs)+1); % this is questionable
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
        u_i = f * dot(disttolmk_k_N',ihatC_k_N) / dot(disttolmk_k_N',khatC_k_N) + u0;
        v_i = f * dot(disttolmk_k_N',jhatC_k_N) / dot(disttolmk_k_N',khatC_k_N) + v0;

       % Check whether landmark is in camera FOV
        if u_i >= 0 && u_i <= umax && dot(disttolmk_k_N',khatC_k_N) > 0
           % Check whether landmark is facing satellite
            if dot(lmkipos_k_N',khatC_k_N) < 0
               % If so, append [timestamp  lmk_id  u   v] to y_table_ideal
                y_table_ideal = [y_table_ideal;...
                                 t_k i u_i v_i]; %#ok<*AGROW> 
            end
        end

    end

end


% Position state history
posStateHist = figure;
subplot(3,1,1);
plot(t,X(:,1),'Color',mlc(1),'DisplayName','$x$');
labels(gca,{'Time [k]','x [km]'},'Position State History');
subplot(3,1,2);
plot(t,X(:,2),'Color',mlc(2),'DisplayName','$y$');
labels(gca,{'Time [k]','y [km]'},'');
subplot(3,1,3);
plot(t,X(:,3),'Color',mlc(3),'DisplayName','$z$');
labels(gca,{'Time [k]','z [km]'},'');
fixfig(posStateHist);
% Velocity state history
velStateHist = figure;
subplot(3,1,1);
plot(t,X(:,4),'Color',mlc(1),'DisplayName','$\dot{x}$');
labels(gca,{'Time [k]','$\mathrm{\dot{x}}$ [km/s]'},'Velocity State History');
subplot(3,1,2);
plot(t,X(:,5),'Color',mlc(2),'DisplayName','$\dot{y}$');
labels(gca,{'Time [k]','$\mathrm{\dot{y}}$ [km/s]'},'');
subplot(3,1,3);
plot(t,X(:,6),'Color',mlc(3),'DisplayName','$\dot{z}$');
labels(gca,{'Time [k]','$\mathrm{\dot{z}}$ [km/s]'},'');
fixfig(velStateHist);


%% Plot Bennu with landmarks
if ~ANIMATE_FLAG
    return
end

f = figure;

% == Plot Bennu ==
xlmks = pos_lmks_A(1,:);
ylmks = pos_lmks_A(2,:);
zlmks = pos_lmks_A(3,:);
[k1,av1] = convhull(xlmks,ylmks,zlmks); % creates convex connectivity matrix for Bennu's surfaces

bennusurf = trisurf(k1,xlmks,ylmks,zlmks,'FaceColor',[0.4 0.4 0.4]); hold on;
set(gca,'Color',[0 0 0]);
axis vis3d; ax = gca;
material dull; light('Position',[1.5e8,0,0]);

% == Plot landmarks ==
lmks = plot3(pos_lmks_A(1,:),pos_lmks_A(2,:),pos_lmks_A(3,:),'c.','MarkerSize',18);
% == Plot satellite track ==
satTrack = plot3(r0(1),r0(2),r0(3),'-','Color',[0.6 0.6 0.6],'LineWidth',2);
% == Plot LOSbeam ==
LOSbeam = plot3([0 r0(1)],[0 r0(2)],[0 r0(3)],'-g','LineWidth',2);
% == Plot satellite ==
sat = plot3(r0(1),r0(2),r0(3),'.','Color',mlc(3),'MarkerSize',30);

set(ax,'XLim',[-1.5 1.5],'YLim',[-1.5 1.5],'ZLim',[-1.5 1.5]);

% Animate
dtheta = dt_int*(2*pi)/(T_A*3600); C = rotZ(dtheta);
for k=1:length(tspan)
    % Rotate Bennu
    bennusurf.Vertices = (C*bennusurf.Vertices')';
    % Rotate landmarks
    lmkdata = C*[lmks.XData; lmks.YData; lmks.ZData];
    lmks.XData = lmkdata(1,:); lmks.YData = lmkdata(2,:);
    % Propagate tracker (orbit line)
    set(satTrack,'XData',[satTrack.XData X(k,1)],'YData',[satTrack.YData X(k,2)],'ZData',[satTrack.ZData X(k,3)]);
    % Propagate satellite
    set(sat,'XData',X(k,1),'YData',X(k,2),'ZData',X(k,3));
    % Propagate LOSbeam
    set(LOSbeam,'XData',[0 X(k,1)],'YData',[0 X(k,2)],'ZData',[0 X(k,3)]);

    pause(0.1);
end
