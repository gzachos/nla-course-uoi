function [X] = solve_system(A, B, n)
	% Perform Cholesky row-decomposition
	L = cholesky(A, n);
	L_trn = L';
	% A = L * L_trn;  % Verify decomposition
	% L * Y = B, solve for Y using forward substitution
	Y = forward_substitution(L, B, n);
	% L^T * X = Y, solve for X using back substitution
	X = back_substitution(L_trn, Y, n);
endfunction
