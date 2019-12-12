function [L] = cholesky(A,n)
	L = zeros(n,n,"single");

	for i = 1:n
		for j = 1:i-1
			tmp_sum = A(i,j);
			for k = 1:j-1
				tmp_sum -= L(i,k) * L(j,k);
			endfor
			L(i,j) = tmp_sum / L(j,j);
		endfor
		tmp_sum = A(i,i);
		for j = 1:i-1
			tmp_sum -= L(i,j) * L(i,j);
		endfor
		L(i,i) = sqrt(tmp_sum);
	endfor
endfunction
