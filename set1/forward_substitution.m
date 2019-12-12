function [y] = forward_substitution(l,b,n)
	y = zeros(10,1,"single");

	for i = 1:n
		tmp_sum = 0.0;
		for k = 1:i-1
			tmp_sum += l(i,k) * y(k);
		endfor
		y(i) = (b(i) - tmp_sum) / l(i,i);
	endfor
endfunction
