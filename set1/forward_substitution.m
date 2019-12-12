function [Y] = forward_substitution(L,B,n)
	Y = zeros(10,1,"single");

	for i = 1:n
		tmp_sum = 0.0;
		t = 1;
		if (i >= 3)
			t = i-2;
		endif
		for k = t:i-1
			tmp_sum += L(i,k) * Y(k);
		endfor
		Y(i) = (B(i) - tmp_sum) / L(i,i);
	endfor
endfunction
