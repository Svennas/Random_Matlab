f = @(x) (x+1).^2.*cos((2*x+1)./(x-3.3));
a = 


function [In,t] = adaptiveSimpson(f,a,b,tol)
%%
%  Adaptive Simpson's rule to approximate a definite integral. This uses
%  a posterioi error analysis to recursively apply Simpson's rule where it
%  is needed
%
%  INPUT:
%
%     f   - integrand function
%     a   - lower bound of interval
%     b   - upper bound of the interval
%     tol - error tolerance for adaptivity
%
%  OUTPUT:
%
%     In - adaptive approximation of the definite integral to
%          within the error tolerance
%     t  - vector of the adapted quadrature nodes
%
%%
%%
%  use recursive bisection with error estimation to compute the integral approximation
   m      = 0.5*(b+a); % find the current midpoint
   [In,t] = do_integral(a,f(a),b,f(b),m,f(m),tol);
%%
%  this is a recursive function within MATLAB so it calls itself
      function [In,t] = do_integral(a,fa,b,fb,m,fm,tol)
      %%
      %  need the two midpoints of the sub-intervals and the function
      %  evalutions for the recursion
         xL  = 0.5*(a+m);
         fxL = f(xL);
         xR = 0.5*(b+m);
         fxR = f(xR);
      %%
      %  save the five nodes at the current level of the recursion
         t = [a;xL;m;xR;b];
      %%
      %  get the h value for the whole interval
         h = 0.5*(b-a);
      %%
      %  compute the Simpson's rule on the whole interval
         S_coarse = h/3*(fa + 4*fm + fb);
      %%
      %  compute the Simpson's rule on the two subintervals
         S_Left  = h/6*(fa + 4*fxL + fm);
         S_Right = h/6*(fm + 4*fxR + fb);
         S_fine  = S_Left + S_Right;
      %%
      %  error estimate
         E = S_coarse - S_fine;
      %%
      %  check against the user set tolerance
         if abs(E) < 10*tol
      % exit when tolerence is met
            In = S_fine;
         else
      %%
      %  bisect again if the tolernce is note met
            [IL,tL] = do_integral(a,fa,m,fm,xL,fxL,tol);
            [IR,tR] = do_integral(m,fm,b,fb,xR,fxR,tol);
            In = IL + IR;
      % append the node edges together (this avoids duplicates)
            t = [tL;tR(2:end)];
         end
      end
end
