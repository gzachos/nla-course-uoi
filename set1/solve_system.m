function [X] = solve_system(A, B, n)
        L = cholesky(A, n);
	% LR = chol(A, "lower");
	% L-LR
        L_trn = L';
        A = L * L_trn;
        Y = forward_substitution(L, B, n);
        X = back_substitution(L_trn, Y, n);
endfunction
