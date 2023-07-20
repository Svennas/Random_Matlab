format short;

C = [4 -2 1; 3 2 -6; 1 -1 1];
B = [2 0 4 3; -2 0 2 -13; 1 15 2 -4.5; -4 5 -7 10];

disp(lufact(B));

function [L,U,P] = lufact(A)
%%
%  Computes the LU factorization of the matrix A. The factorization is done
%  in-place and then the L and U matrices are extracted at the end for
%  output. The factorization is PA = LU
%
%  INPUT: 
%    A - n by n square, non-singular matrix
%
%  OUTPUT:
%    L - lower triangular matrix with ones along the main diagonal
%    U - upper triangular matrix
%    P - permutation matrix for pivoting
%%
%  get the size of the system
   n  = length(A);
%%
%  initialize the pivoting matrix
   P = eye(n);
%%
%  Vectorized in-place LU factorization (with row pivoting) that keeps
%  track of the total permutations by scrambleing the matrix P

   for j = 1:n-1
%%
% Find the index of the largest pivot element in the current column
      
      v = A(:,j);

      [maxVal, k] = max(v(j : n, :));
         
      k = k + j - 1; 
  
%%
% Swap the rows within the in-place array as well as the permutation matrix P
      A([j k],:) = A([k j],:);
      P([j k],:) = P([k j],:);
%%
% Perform the in-place elimination and save the new column of L
      i = j+1:n; % indices for the "active" matrix portion
      A(i,j) = A(i,j)/A(j,j); 
      A(i,i) = A(i,i) - A(i,j)*A(j,i); 
     % A, return
     
   end
%%
% Extract L and U from the in-place form
   U = triu(A);
   L = eye(n) + tril(A,-1);
   disp("U");
   disp(U);
   disp("L");
   disp(L);
   %disp("A");
   %disp(A);
   disp("P");
   disp(P);
end
