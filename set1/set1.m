#!/usr/bin/env octave

n = 1000

A1 = zeros(n,n,"single");
A2 = zeros(n,n,"single");
B1 = zeros(n,1,"single");
B2 =  ones(n,1,"single");

% Initialize B1
B1(1,1) = B1(n,1) = 3;
B1(2,1) = B1(n-1,1) = -1;

% Initialize B2
B2(1,1) = B2(n,1) = 4;
B2(2,1) = B2(n-1,1) = 0;

% Initialize A1, A2
for i = 1:n
        for j = 1:n
                if (i == j)
                        A1(i,i) = 6;
                        A2(i,i) = 7;
                elseif (i+j == (i*2)-1 || i+j == (j*2)-1)
                        A1(i,j) = A2(i,j) = -4;
                elseif (i+j == (i*2)-2 || i+j == (j*2)-2)
                        A1(i,j) = A2(i,j) = 1;
                endif
        endfor
endfor


X1 = solve_system(A1, B1, n)

X2 = solve_system(A2, B2, n)

filename1 = sprintf('x1_%d.mat', n);
filename2 = sprintf('x2_%d.mat', n);

save_1d_matrix(filename1, X1, n);
save_1d_matrix(filename2, X2, n);

% save file1 X1
% save file2 X2

