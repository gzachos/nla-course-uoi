#!/usr/bin/env octave

n = 10

A1 = zeros(n,n);

for i = 1:n
        for j = 1:n
                if (i == j)
                        A1(i,i) = 6;
                elseif (i+j == (i*2)-1 || i+j == (j*2)-1)
                        A1(i,j) = -4;
                elseif (i+j == (i*2)-2 || i+j == (j*2)-2)
                        A1(i,j) = 1;
                endif
        endfor
endfor

A1

L = cholesky(A1, n)

