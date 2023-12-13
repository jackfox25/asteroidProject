load('x_noisy_MC40.mat');

ff = figure(1); hold off; clf;
plt0 = plot3(x_nom(1,1),x_nom(1,2),x_nom(1,3),'Color',0.7*mlc(1)); hold on;
plt0h = plot3(x_nom(1,1),x_nom(1,2),x_nom(1,3),'.','Color',mlc(1),'MarkerSize',35);
hold on; axis equal; axis vis3d;

% == Plot Bennu ==
xlmks = pos_lmks_A(1,:);
ylmks = pos_lmks_A(2,:);
zlmks = pos_lmks_A(3,:);
[k1,av1] = convhull(xlmks,ylmks,zlmks); % creates convex connectivity matrix for Bennu's surfaces

bennusurf = trisurf(k1,xlmks,ylmks,zlmks,'FaceColor',[0.4 0.4 0.4]); hold on;
set(gca,'Color',[0 0 0]);
axis vis3d; ax = gca;
material shiny; light('Position',[1.5e8,0,0]);

% == Plot landmarks ==
lmks = plot3(pos_lmks_A(1,:),pos_lmks_A(2,:),pos_lmks_A(3,:),'c.','MarkerSize',18);

set(ax,'XLim',[-1.5 1.5],'YLim',[-1.5 1.5],'ZLim',[-1.5 1.5]);



plt1 = plot3(x_noisy_MC(1,1,1),x_noisy_MC(1,2,1),x_noisy_MC(1,3,1),'Color',0.7*mlc(2));
plt2 = plot3(x_noisy_MC(1,1,2),x_noisy_MC(1,2,2),x_noisy_MC(1,3,2),'Color',0.7*mlc(3));
plt3 = plot3(x_noisy_MC(1,1,3),x_noisy_MC(1,2,3),x_noisy_MC(1,3,3),'Color',0.7*mlc(4));

plt1h = plot3(x_noisy_MC(1,1,1),x_noisy_MC(1,2,1),x_noisy_MC(1,3,1),'.','Color',mlc(2),'MarkerSize',30);
plt2h = plot3(x_noisy_MC(1,1,2),x_noisy_MC(1,2,2),x_noisy_MC(1,3,2),'.','Color',mlc(3),'MarkerSize',30);
plt3h = plot3(x_noisy_MC(1,1,3),x_noisy_MC(1,2,3),x_noisy_MC(1,3,3),'.','Color',mlc(4),'MarkerSize',30);

view(-90,0);
labels(gca,{'$X$ [km]','$Y$ [km]','$Z$ [km]'},'');
camzoom(1.5); fixfig(ff);

dtheta = dt_int*(2*pi)/(T_A); C = rotZ(dtheta);

ll = legend('NOMINAL','interpreter','latex');
ll.Color = [1 1 1];

for i=2:25:size(x_nom,1)

    set(plt0h,'XData',x_nom(i,1),'YData',x_nom(i,2),'ZData',x_nom(i,3));
    set(plt1h,'XData',x_noisy_MC(i,1,1),'YData',x_noisy_MC(i,2,1),'ZData',x_noisy_MC(i,3,1));
    set(plt2h,'XData',x_noisy_MC(i,1,2),'YData',x_noisy_MC(i,2,2),'ZData',x_noisy_MC(i,3,2));
    set(plt3h,'XData',x_noisy_MC(i,1,3),'YData',x_noisy_MC(i,2,3),'ZData',x_noisy_MC(i,3,3));

    set(plt0,'XData',[plt0.XData x_nom(i,1)],'YData',[plt0.YData x_nom(i,2)],'ZData',[plt0.ZData x_nom(i,3)]);
    set(plt1,'XData',[plt1.XData x_noisy_MC(i,1,1)],'YData',[plt1.YData x_noisy_MC(i,2,1)],'ZData',[plt1.ZData x_noisy_MC(i,3,1)]);
    set(plt2,'XData',[plt2.XData x_noisy_MC(i,1,2)],'YData',[plt2.YData x_noisy_MC(i,2,2)],'ZData',[plt2.ZData x_noisy_MC(i,3,2)]);
    set(plt3,'XData',[plt3.XData x_noisy_MC(i,1,3)],'YData',[plt3.YData x_noisy_MC(i,2,3)],'ZData',[plt3.ZData x_noisy_MC(i,3,3)]);

    
    % Rotate Bennu
    bennusurf.Vertices = (C*bennusurf.Vertices')';
    % Rotate landmarks
    lmkdata = C*[lmks.XData; lmks.YData; lmks.ZData];
    lmks.XData = lmkdata(1,:); lmks.YData = lmkdata(2,:);

    pause(0.1);

end