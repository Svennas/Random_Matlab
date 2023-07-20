function plot_composite_quad(f,t)
%%
%  Script that will plot a function and the interval boundaries for use
%  with a compositie quadrature rule
%
%  INPUT:
%
%     f - integrand function
%     t - set of nodes used for the composite quadrature rule
%
%  OUTPUT:
%
%     NONE - Produces plot to screen
%
%%
   clf
   p1 = fplot(f,[t(1) t(end)],'-k','LineWidth',2);
   hold on
   p2 = stem(t,f(t),'fill','color',[0.9100    0.4100    0.1700],'MarkerSize',9,'LineWidth',1.5);
   set(gcf,'Position',[500, 60, 1250, 1250])
   set(gca,'FontSize',16);
   legend([p1 p2],{'$(x+1)^2\cos\left(\frac{2x+1}{x-3.3}\right)$','Quadrature nodes'}...
          ,'interpreter','latex','fontsize',24,'Location','northwest')
end