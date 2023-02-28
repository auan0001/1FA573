clear
% Samples and walks
N = 10;

M = 200;
Nx = 1:10

for samples = 1:Nx(i)
  % x = zeros(1,Nx(i)+1);
  % y = zeros(1,Nx(i)+1);
  x(1) = 0;
  y(1) = 0;
  for step = 1:M+1 
    if sign(randn) > 0
      x(step+1) = x(step) + sign(randn);
    else
      y(step+1) = y(step) + sign(randn);
    end
    % plot(x, y, 'b');
    % grid on
    % hold on
  end
  R2(samples) = sqrt(sum(x.^2+y.^2)^2)
end
  

% hold off
% subplot(2,1,1)
Nx
plot(Nx,sort(R2/N), 'r.')
% title('Distance vs sample size')
% grid on
% xlabel('N^{1/2}')
% ylabel('<R^2>^{1/2}')
% subplot(2,1,2)
% loglog((mean(R2,2)), (N_incr))
% title('Loglog distance vs sample size')
% xlabel('N^{1/2}')
% ylabel('<R^2>^{1/2}')
% grid on
