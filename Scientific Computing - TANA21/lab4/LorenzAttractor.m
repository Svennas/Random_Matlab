T = 50;



function LorenzAttractor(T)
%%
%  Nonlinear ODE example from nonlinear dynamics to solve the
%  Lorenz system using a Runge-Kutta 4 time integrator
%
%      d(y1)/dt = 10*(y2-y1)
%      d(y2)/dt = 28*y1 - y2 - y1*y3
%      d(y3)/dt = y1*y2 - 8/3*By3
%
%  INPUT:
%
%     T - final time
%
%  OUTPUT:
%
%     NONE - this makes a movie
%
%%
%  adhoc tactic to set the time step
   num_steps = 2500;
   dt = T/num_steps;
%%
% set the initial condition and allocate memory for the solution vector
   y   = zeros(3,num_steps);
   ini = 5*ones(3,1);
%%
% save the initial condition and enter the Runge-Kutta loop to solve the ODE
   y(:,1) = ini;
   for i = 2:num_steps
      k1 = RHS(y(:,i-1));
      k2 = RHS(y(:,i-1) + 0.5*dt.*k1);
      k3 = RHS(y(:,i-1) + 0.5*dt.*k2);
      k4 = RHS(y(:,i-1) + dt.*k3);
      y(:,i) = y(:,i-1) + (dt/6).*(k1 + 2.*k2 + 2.*k3 + k4);
   end
%%
% visualize the solution progessively to make a movie
   clf;
   hold on;
   xlim([-20 20]);
   ylim([-30 30]);
   zlim([0 50]);
   grid on;
   for n = 1:(length(y)/7)
      plot3(y(1,1:7*n),y(2,1:7*n),y(3,1:7*n),'b');
      drawnow
      view(-37.5-n, 24)
      pause(0.02)
   end
end

%%%
%  auxiliary function to get the right hand side
%%%
 function dy = RHS(y)
%%
% evaluate the right hand side of the nonlinear ODE for the Lorenz
% attractor
   dy    = zeros(3,1);
   dy(1) = 10*(y(2) - y(1));
   dy(2) = 28*y(1) - y(2) - y(1)*y(3);
   dy(3) = y(1)*y(2) - (8/3)*y(3);
 end