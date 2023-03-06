clear

N = 8;
print = 1;

for i = 1:N
  for j = 1:N
    A(i,j) = (j-1)+(i-1)*(N)+1;
  end
end

% Spin matrix
A = (-1).^randi(2,N,N);

heatmap(A);

for i = 1:N
  for j = 1:N
    E(i,j) = 2*A(i,j) ...
      *(A(i, mod(j-2, N)+1) ...
      + A(mod(i-2, N)+1, j) ...
      + A(i ,mod(j, N)+1) ...
      + A(mod(i, N)+1, j));
  end
end


