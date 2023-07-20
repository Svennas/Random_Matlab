%z = 1;%@(k) k;
fk = @(x) (x.^3) + 2*x - 3;
%fk = @(x) x.^2.;
n = 3;
%xVals = [-1 0 1 1.5];
%yVals = [-6 -3 0 3.375];
%xVals = [1 2 3 4 5 6];
%yVals = [1 1 2 6 24 120];
xVals = [0.0000, 0.0016, 0.0072, 0.0185, 0.0279, 0.0462, 0.0625, 0.0919, 0.1109, 0.1354, 0.2500, 0.3600, 0.4900, 0.6400, 0.8100, 1.0000];
yVals = [0.0000, 0.0070, 0.0146, 0.0228, 0.0275, 0.0344, 0.0390, 0.0454, 0.0486, 0.0519, 0.0594, 0.0592, 0.0535, 0.0420, 0.0246, 0.0000];
M = 150;




%disp(Lagrange_interpolation(6, xVals, yVals));
%disp(gamma(6));

build_interpolation(M, xVals, yVals);

function [X, L] = build_interpolation(M, xVals, yVals)
    X = linspace(0, 1, 150);
    
    for k = 1:M
        L(k) = Lagrange_interpolation(X(k), xVals, yVals);
    end 
   
    s = spline(xVals,yVals,X);
    
    plot(X, s,'bo','MarkerFaceColor','b','MarkerSize',8);
    %disp(norm(fk(X) - L, inf));
    %disp((X-1)/2);
end

function [pnx] = Lagrange_interpolation(z, xVals, yVals)
%%
%  Implementation to evaluate an interpolating polynomial p_n(x) 
%  at the point x = z. The polynomial uses the standard Lagrange
%  basis functions.
%
%  INPUT:
%
%     z      - 1x1 value to evaluate
%     fk     - A function
%     xVals  - Matrix with all x values
%  OUTPUT:
%
%     pnx - value of the polynomial interpolant at x = z
%

%%
%  compute the polynomial interpolation sum evaluated at x = z
   n = length(xVals);
    
   x = xVals;
   y = yVals;
   l = ones(1, n);
   pnx = 0;
   
   for i = 1:n
      for j = 1:n
          if x(j) ~= x(i)
              prod = (z - x(j)) / (x(i) - x(j));
              l(i) = l(i) * prod;
          end
      end
      pnx = pnx + (y(i) * l(i));
   end
end

