i = 1;

x = [-1 0 1];

xi = x(i);

disp(x(i));

%prod = sumprod((1 - x(j)) / (xi - x(j)), j, 1, 2);

%disp(prod);

product = 1;

n = 10;

for j = 1:n
   prod = (x - x(j)) / (xi - x(j));
   product = product * prod;
end

disp(product);