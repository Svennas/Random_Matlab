z = 0;
fx = 
intmin= 0;
intmax= 1;

Lagrange_interpolation(z, fx, n, intMin, intMax)

function [pnx] = Lagrange_interpolation(z, fx, n, intMin, intMax)
%%
%  Implementation to evaluate an interpolating polynomial p_n(x) 
%  at the point x = z. The polynomial uses the standard Lagrange
%  basis functions.
%
%  INPUT:
%
%     z    - 1x1 value to evaluate
%     ???
%
%  OUTPUT:
%
%     pnx - value of the polynomial interpolant at x = z
%

%%
%  compute the polynomial interpolation sum evaluated at x = z
   x = zeros(1, n);
   y = zeros(1, n);
   l = ones(1, n);
   step = (abs(intMin) + abs(intMax)) / n;
   pnx = 0;
   
   for i = 1:n
       x(i) = intMin + (step * i);
       y(i) = fx(x(i));
      for j = 1:n
          if j ~= i
              prod = (z - x(j)) / (x(i) - x(j));
          end
          l(i) = l(i) * prod;
         
      end

      pnx = pnx + (y(i) * l(i));
   end
end
