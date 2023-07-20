%format long e

n = 20;
f = @(x) (x+1).^2.*cos((2*x+1)./(x-3.3));
a = 0;
b = 3;


I = integral(f,0,3);
disp("I");
disp(I);

[t, x20] = SimpsonsRule(f, a, b, n, I);

disp(t);
%SimpsonsRule(f, a, b, n, I);



function [S,x] = SimpsonsRule(f,a,b,n, I)
%%
%  Implementation of the (possibly) composite Simpsons's quadrature
%  rule to approximate the definite integral
%
%           b
%          /\
%          |
%          | f(x) dx
%          |
%         \/
%         a
%
%  INPUT:
%
%     f    - anonymous function integrand
%     a    - lower integration bound
%     b    - upper integration bound
%     n    - number of intervals (must be an even integer)
%
%  OUTPUT:
%
%     In - approximation to the definite integral on the given interval [a,b]
%     x  - set of n+1 quadrature nodes
%
%%
%  Simpson's rule needs at least two intervals
   if n == 1
      error('Simpsons rule requires n >= 2');
   end
%%
%  check if the number of intervals is even
   if mod(n,2) ~= 0
      error('The number of intervals n must be even');
   end
%%
%  Find the interval spacing h
   h = (b - a) / n;
%%
%  build the set of uniformly spaced quadrature nodes
   x = zeros(n+1,1);

%%
          
%  Initialize the integral value to zero
   S = 0.0;
   sum = 0;
%%
%  Approximate the integral with (composite) Simpson's rule
   for i = 1:n+1
       if i == 1
           x(i) = a;
           sum = f(x(i));
       elseif i == n+1
           x(i) = b;
           sum = sum + f(x(i));
       else
           x(i) = a + (i-1)*h;
           if mod(i, 2) == 0
               sum = sum + (4*f(x(i)));
           else
               sum = sum + (2*f(x(i)));
           end
       end
       
   end
   S = (h/3)*sum;
   disp("S");
   disp(S);
   disp("E");
   disp(abs(I - S));
   
   plot_composite_quad(f,x)
   
end

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