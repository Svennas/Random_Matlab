function[func, absdiff, eoc] = analyze(n, m)

    format short e;
    
    func = 3 * n - cos(2*pi*n);
    
    absdiff = abs((1/6) - n);
    
    eoc = (log(absdiff))/log(abs((1/6) - m));
    
    
    disp(func);
    
    disp(absdiff);
    
    disp(eoc);

end