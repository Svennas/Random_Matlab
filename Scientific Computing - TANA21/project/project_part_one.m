A = full(gallery('tridiag', 10, -1, 2, -1));
b = [2.5; 0; 0; 0; 0; 0; 0; 0; 0; 2.5];
x0 = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
tol = 10e-09;

% Vector that saves the amount of iterations
ints = zeros(1, 1);
% Vector that saves the tested weights
weights = zeros(1, 1);
i = 1;
w = 1.05;

% Create vectors with all weights and all the amount of iterations
while w < 2
   ints(i) = SOR(A, b, x0, tol, w);
   weights(i) = w;
   
   i = i + 1;
   w = w + 0.05;
end

% Plot the relation between weight and iterations
plot(weights, ints, 'LineWidth', 2, 'Marker', 'd');
grid on;
title('Weight to iterations relation')
xlabel('Weight')
ylabel('Iterations')

function [x] = SOR(A,b,x0,tol,w)
%%
%  The Gauss-Seidel iterative method known as
%  Successive Over-relaxation (SOR), which approximates the solution of a
%  linear system Ax=b up to a user defined tolerance
%
%  INPUT: 
%    A   - n by n square, non-singular matrix
%    b   - n by 1 right hand side vector
%    x0  - n by 1 vector containing that initial guess for the iteration
%    tol - user set tolerance for the stopping condition in the iteration 
%    w   - user set weight 
%  OUTPUT:
%    x - n by 1 vector containing the iterative solution
%    k - number of iterations
%%
%  get the system size
   n = length(A);
%%
%  save the initial guess
   x = x0;
   %disp(length(x));
   %disp(length(A));
%%
%  Gauss-Seidel iteration which overwrites the current approximate solution
%  with the new approximate solution

   max = 500;
   for k = 1:max
       for i = 1:n
           sum = 0;
           for j = 1:n
               if j ~= i
                   sum = sum + (A(i,j) * x(j));
               end
               
           end
           x(i) = (((b(i) - sum)/A(i,i)) * w) + (1 - w)*x(i);
       end
        
           r = b - (A*x);
           if norm(r, 2) < tol * norm(b, 2)
               disp(k);
               break
           end
   end
end