function airfoil_image
%%
%  This draws the parametric curves for a NACA0012 airfoil
%
   t = linspace(0,1,150);
% sample the top parametric curve
   topx = (1 - t).^2;
   topy = 0.6*(0.2969*(1-t) - 0.1260*(1-t).^2 - 0.3516*(1-t).^4 +0.2843*(1-t).^6-0.1036*(1-t).^8);
% sample the bottom parametric curve
   botx = t.^2;
   boty = -0.6*(0.2969*t - 0.1260*t.^2 - 0.3516*t.^4 +0.2843*t.^6-0.1036*t.^8);
% plot the top and bottom curves
   plot(topx,topy,'-k',botx,boty,'-k','LineWidth',1.5)
   hold on
% set limits and size of figure
   set(gcf,'Position',[100, 80, 1850, 500])
   set(gca,'FontSize',16);
   xlim([-0.025 1.025])
   ylim([-0.1 0.1])
% fill with gray color for better visualization
   h1 = fill(topx,topy,'k','LineStyle','none');
   set(h1,'facealpha',0.05);
   h2 = fill(botx,boty,'k','LineStyle','none');
   set(h2,'facealpha',0.05);
%
end