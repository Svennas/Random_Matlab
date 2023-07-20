C = full(gallery('tridiag', 10, -1, 2, -1));
b = [2.5; 0; 0; 0; 0; 0; 0; 0; 0; 2.5];
x0 = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
tol = 10e-09;

disp(gSeidel(C, b, x0, tol));


function [x] = gSeidel(A,b,x0,tol)
%%
%  Gauss-Seidel iterative method to approximate the solution of a
%  linear system Ax=b up to a user defined tolerance
%
%  INPUT: 
%    A   - n by n square, non-singular matrix
%    b   - n by 1 right hand side vector
%    x0  - n by 1 vector containing that initial guess for the iteration
%    tol - user set tolerance for the stopping condition in the iteration 
%
%  OUTPUT:
%    x - n by 1 vector containing the iterative solution
%    k - number of iterations
%%
%  get the system size
   n = length(A);
%%
%  save the initial guess
   x = x0;
%%
%  Gauss-Seidel iteration which overwrites the current approximate solution
%  with the new approximate solution

   max = 250;
   for k = 1:max
       xp = x;
       for i = 1:n
           sum = 0;
           for j = 1:n
               if j ~= i
                   sum = sum + (A(i,j) * x(j));
               end
               x(i) = (b(i) - sum)/A(i,i);
           end
       end
           r = (x - xp);
           if ((norm(r, 2) / norm(x, 2)) < tol)
               disp(k);
               break
           end
   end
%
end



