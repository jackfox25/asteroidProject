%% PROBLEM 1 %%
% Noiseless nonlinear simulation

clear; clc; close all;
ANIMATE_FLAG = false;

load('orbitdetermination-finalproj_data_2023_11_14.mat');
simParameters; % load in constants

tspan = t0:dt_int:tf;
tol=1.e-14;
OPTIONS = odeset('RelTol',3*tol,'AbsTol',tol);

x0 = [r0; rdot0];
processNoise = zeros(3,1);
[t,X] = ode45(@satDynamicModel,tspan,x0,OPTIONS,processNoise);

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
