function [y] = back_substitution(l,b,n)
	y = zeros(10,1,"single");

	for i = n:-1:1
		tmp_sum = 0.0;
		for k = n:-1:i+1
			tmp_sum += l(i,k) * y(k);
		endfor
		y(i) = (b(i) - tmp_sum) / l(i,i);
	endfor
endfunction
