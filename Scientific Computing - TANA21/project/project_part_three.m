%%% Create the function for the exact solution
c2 = (8 - 12*(sin(log(2))) - 4*(cos(log(2)))) *( 1/70);
c1 = (11/10) -c2 ;
y = @(x) c1 * x + c2/(x*x) - (3/10)*sin(log(x)) - (1/10)*cos(log(x));
%%%
% Define the tolerance and weight for SOR 
tol = 10e-09;
w = 1.58;

%%%%%
% The sets of nodes (should be changed during testing)
n = 80;
% S is the vector for our solution
S = polli( n, tol, w);
% Defines h so the exact solution can be calculated
h = 1/(n-1);
% Creates a vector with the exact solution
exact = zeros(n, 1);
    for i = 0:n-1
        exact(i+1) = y(1+(i*h)); 
    end
% ints (iterations) is a vector with all h used 
ints = (1:h:2)';
% Plots the exact solution (e_sol) and our estimated solution (M)
plot(S, ints, exact, ints, 'LineWidth', 2);
grid on;
title('Approximated function (red line) and exact function (blue line)')
xlabel('x')
xlabel('x')
ylabel('Function values')
%%%%%



function [S] = polli(n, tol, w)
%%% inputs
% n - amount of nodes
% tol - the tolerance for SOR
% w - the weight for SOR
%%% outputs
% S - the solution vector after using SOR
%%%
    % Define h
    h = 1/(n-1);
    
    % Create function to use for vector b
    fb = @(x) -(sin(log(x))) / 2;
  
    % Create vector with all values for h
    int = (1:h:2)';
    
    % Makes P1
    vP1 = zeros(n, 1);
    fP1 = @(x) 2/x;
    for i = 1:n
        vP1(i) = fP1(int(i));
    end
    P1 = diag(vP1);
    
    % Makes P2
    fk = @(x) 2/(x.^2);
    vP2 = zeros(n, 1);
    for i = 1:n
        vP2(i) = fk(int(i));
    end
    P2 = diag(vP2);
    
     % Makes D1
    D1 = full(gallery('tridiag', n, -0.5, 0, 0.5));
    D1(1, 1) = -1;
    D1(1, 2) = -1;
    D1(n, n) = 1;
    D1(n, n-1) = -1;
    D1 = D1/h;
    
    % Makes D2
    D2 = full(gallery('tridiag', n, 1, -2, 1));
    
    D2(1, 1) = 2;
    D2(1, 2) = -5;
    D2(1, 3) = 4;
    D2(1, 4) = -1;
    
    D2(n, n) = 2;
    D2(n, n-1) = -5;
    D2(n, n-2) = 4;
    D2(n, n-3) = -1;
    
    D2 = D2/(h*h);
    
    % Makes b (from Ax=b)
    b = zeros(n, 1);
    for i = 0:n-1
       b(i+1) = fb(1+i*h); 
    end
    b(1) = 1;
    b(n) = 2;
    
    % Makes A (from Ax=b)
    A = (D2 + (P1*D1) -P2);
    A(1, :) = 0;
    A(1, 1) = 1;
    A(n, :) = 0;
    A(n, n) = 1;
    
    % Create vector with initial guess
    x0 = ones(n, 1);
    
    % Create our solution S
    S = SOR(A, b, x0, tol, w);
end

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
   %disp(length(x));
   %disp(length(A));
%%
%  Gauss-Seidel iteration which overwrites the current approximate solution
%  with the new approximate solution

   max = 250;
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
%
disp("This is our solution")
disp(x);
end