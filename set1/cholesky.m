
function [l] = cholesky(a,n)
	l = zeros(n,n);

	for i = 1:n
		for j = 1:i-1
			tmp_sum = a(i,j);
			for k = 1:j-1
				tmp_sum -= l(i,k) * l(j,k);
			endfor
			l(i,j) = tmp_sum / l(j,j);
		endfor
		tmp_sum = a(i,i);
		for j = 1:i-1
			tmp_sum -= l(i,j) * l(i,j);
		endfor
		l(i,i) = sqrt(tmp_sum);
	endfor
endfunction

