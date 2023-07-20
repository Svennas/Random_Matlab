func = @(x) 3*x - cos(2*pi*x);
x1 = 0.4;
toll = 1e-09;

disp(steffe(func, x1, toll));

function [xVals,iter] = steffe(f,x0,tol)
%%
%  Implementation of the Newton-Raphson method to approximate 
%  the root of a nonlinear function f(x)
%
%  INPUT: 
%    f      - a function f(x)
%    fprime - first derivative of the function f(x)
%    x0     - initial guess of the root location
%    tol    - error tolerance to stop the iteration
%
%  OUTPUT:
%    xVals  - sequence of approximate values for the root of the function f(x)
%    iter   - number of iterations it took to achieve user set tolerance
%%
%  initialize a vector to store the sequence of guesses
   xVals = zeros(16,1);
%%
%  save the initial guess
   xVals(1) = x0;
%%
%  perform Steffensen's for a fixed number of iterations
   for k = 1:15
%%
%  update to the next value
      fk  = f(xVals(k));
      xVals(k+1) = xVals(k) - (fk*fk) / (f(xVals(k)+fk) - fk);
%%
% check the kill condition
      stopCond = abs(f(xVals(k+1)));
      if stopCond < tol
         xVals = xVals(1:k+1);
         iter  = k+1;
         break;
      end
   end
   disp("Iterations");
   disp(iter);
end