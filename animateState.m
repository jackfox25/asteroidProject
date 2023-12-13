load('x_noisy_MC40.mat');

ff = figure(1); hold off; clf;
plt0 = plot3(x_nom(1,1),x_nom(1,2),x_nom(1,3),'Color',mlc(1),'LineWidth',3);
hold on; axis equal; axis vis3d;
axis([-0.5 0.5 -1.5 1.5 -1.5 1.5]);


plt1 = plot3(x_noisy_MC(1,1,1),x_noisy_MC(1,2,1),x_noisy_MC(1,3,1),'Color',mlc(2));
plt2 = plot3(x_noisy_MC(1,1,2),x_noisy_MC(1,2,2),x_noisy_MC(1,3,2),'Color',mlc(3));
plt3 = plot3(x_noisy_MC(1,1,3),x_noisy_MC(1,2,3),x_noisy_MC(1,3,3),'Color',mlc(4));

movefig(ff,'full')

for i=2:25:size(x_nom,1)

    set(plt0,'XData',[plt0.XData x_nom(i,1)],'YData',[plt0.YData x_nom(i,2)],'ZData',[plt0.ZData x_nom(i,3)])
    set(plt1,'XData',[plt1.XData x_noisy_MC(i,1,1)],'YData',[plt1.YData x_noisy_MC(i,2,1)],'ZData',[plt1.ZData x_noisy_MC(i,3,1)]);
    set(plt2,'XData',[plt2.XData x_noisy_MC(i,1,2)],'YData',[plt2.YData x_noisy_MC(i,2,2)],'ZData',[plt2.ZData x_noisy_MC(i,3,2)]);
    set(plt3,'XData',[plt3.XData x_noisy_MC(i,1,40)],'YData',[plt3.YData x_noisy_MC(i,2,40)],'ZData',[plt3.ZData x_noisy_MC(i,3,40)]);

    pause(0.1);

end