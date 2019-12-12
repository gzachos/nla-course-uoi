function [L] = cholesky(A,n)
	L = zeros(n,n,"single");
	for i = 1:n
		t1 = 1;
		if (i > 2)
			t1 = i-2;
		endif
		for j = t1:i-1
			tmp_sum = 0.0;
			t2 = 1;
			if (j > 2)
				t2 = j-2;
			endif
			for k = t2:j-1
				tmp_sum += L(i,k) * L(j,k);
			endfor
			L(i,j) = (A(i,j) - tmp_sum) / L(j,j);
		endfor
		tmp_sum = 0.0;
		t1 = 1;
		if (i > 2)
			t1 = i-2;
		endif
		for j = t1:i-1
			tmp_sum += L(i,j) * L(i,j);
		endfor
		L(i,i) = sqrt(A(i,i) - tmp_sum);
	endfor
endfunction
