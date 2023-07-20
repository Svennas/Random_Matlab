function [xVals,iter] = newtonRaphson(f,fprime,x0,tol)
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
   xVals = zeros(10000,1);
%%
%  save the initial guess
   xVals(1) = x0;
   stopCond = 1;
   k = 1;
%%
%  perform Newton-Raphson for a fixed number of iterations
   while (stopCond > tol)
%%
%  update to the next value
      fk  = f(xVals(k));
      dfk = fprime(xVals(k));
      xVals(k+1) = xVals(k) - fk/dfk;
%%
% check the kill condition
      stopCond = abs(f(xVals(k+1)))
      iter  = k+1;
      k = k + 1;
   end
   xVals = xVals(1:iter);
end