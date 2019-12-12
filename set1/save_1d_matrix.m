function save_1d_matrix(filename, A, n)
	fd = fopen(filename, 'w');
	for i = 1:n
		fprintf(fd, "%.5f\n", A(i));
	endfor
	fclose(fd);
endfunction
