C = gallery('circul', -4:4);
S = magic(7);
Dx = toeplitz([6 -4 1 zeros(1,12)]);

disp(C);
%disp(cond(Dx));
disp(eig(C));
%disp(det(Dx));
disp(S);
disp(Dx);