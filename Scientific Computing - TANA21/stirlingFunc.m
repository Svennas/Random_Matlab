format short e

for i = 2:12
   stirlingFunct(i); 
end

function[fact, stirling, absolute, relative] = stirlingFunct(n)
    e = exp(1);
    
    fact = factorial(n);
    
    stirling = (sqrt(2*pi*n))*((n/e)^n);
    
    absolute = fact - stirling;
    
    relative = absolute / fact;
    
    disp("n");
    disp(n);
    disp("Factorial:");
    disp(fact);
    disp("Stirling");
    disp(stirling);
    disp("Absolute");
    disp(absolute);
    disp("Relative");
    disp(relative);
end


