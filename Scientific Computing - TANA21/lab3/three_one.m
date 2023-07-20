f = @(x)exp(sin(4*x));
df = @(x) 4*(cos(4*x)*(exp(sin(4*x))));
format short;
%f = @(x);
%df = @(x)2*x;
derivative (f, df);

function [maxv] = derivative (fx, df)
    N = 10;
    f = zeros(N,1);
    fprimes = f;
    maxv = 0;
    h = 1 / N;

    for i = 0:N
        x = i/N;
        f(i + 1) = fx(x);
        %disp(df(x));
    end
    
    for i = 1:N
        x = (i-1) / N;
        if i == 1
            fprime = (f(i+1) - f(i)) / h;
            %fprime = (4*(f(i+1)) - 3*(f(i)) - f(i+2))/(2*h);
            fprimes(i) = fprime;
            
            if maxv < (abs(df(x) - fprimes(i)))
                maxv = (abs(df(x) - fprimes(i)));
                %disp(maxv);
            end
            
        elseif i == N
            %fprime = (3*f(i) - 4*f(i-1) + f(i-2)) / (2*h);
            fprime = (f(i) - f(i-1)) / h;
            fprimes(i) = fprime;
            
            if maxv < (abs(df(x) - fprimes(i)))
                maxv = (abs(df(x) - fprimes(i)));
                disp("N")
                disp(maxv);
            end
                
        else
            fprime = (f(i + 1) - f(i - 1))/(2*h);
            fprimes(i) = fprime;
            
            if maxv < (abs(df(x) - fprimes(i)))
                maxv = (abs(df(x) - fprimes(i)));
                %disp(maxv);
            end
            
        end
    end
    disp(fprimes);
    disp (maxv);
end