function [Y] = back_substitution(L,B,n)
	Y = zeros(n,1,"single");

	for i = n:-1:1
		tmp_sum = 0.0;
		t = n;
		if (i+2 <= n)
			t = i+2;
		endif
		for k = t:-1:i+1
			tmp_sum += L(i,k) * Y(k);
		endfor
		Y(i) = (B(i) - tmp_sum) / L(i,i);
	endfor
endfunction
