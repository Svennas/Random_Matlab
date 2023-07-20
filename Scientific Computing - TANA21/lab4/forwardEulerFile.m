f = @(t,y) 5./t - 1./(t.^2) - 5.*t*y.^2;
t0 = 1;
T = 25;
y0 = 1;
h = 0.1;

forwardEuler(f,t0,T,y0,h);

function [t,y] = forwardEuler(f,t0,T,y0,h)
%%
%  Forward Euler time integration technique to solve the
%  initial value problem
%
%     y' = f(t,y),  y(t0) = y0
%
%  INPUT:
%
%     f  - right-hand-side function that can depend on
%          the time t and the function y
%     t0 - intial time
%     T  - final time
%     y0 - initial solution value
%     h  - fixed time step size
%
%  OUTPUT:
%
%     t - vector of time values
%     y - solution vector of approximate values
%
%%
%  discretize the time variable and set the problem size
   t = t0:h:T;
   y = zeros(length(t),1);
%%
%  save the intial values
   y(1) = y0;
   t(1) = t0;
%%
%  use forward Euler with the given timestep h to solve the problem and
%  save the information for later analysis
   for j = 2:length(t)
%  this is an explicit method so we always use values from the current solution
      y(j) = y(j-1) + h*f(t(j-1),y(j-1));
   end
   disp(y(end));
end
