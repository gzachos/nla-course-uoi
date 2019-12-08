function [y] = forward_substitution(l,b,n)
	y = zeros(10,1);
	for i = 1:n
		tmp_sum = b(i);
		for k = 1:i-1
			tmp_sum -= l(i,k) * y(k);
		endfor
		y(i) = tmp_sum / l(i,i);
	endfor
endfunction
