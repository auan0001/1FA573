clear
% Samples and walks
N = 5000;
N_incr = linspace(100,N,N/1000);

% Start from origo
% x = zeros(N+1,1);
% y = zeros(N+1,1);
for i = 1:length(N_incr)
  for samples = 1:N_incr(i)
    x = zeros(N_incr(i)+1,1);
    y = zeros(N_incr(i)+1,1);
    for step = 1:N_incr(i) 
      % sign() is -1 or +1
      x(step+1) = x(step) + sign(randn);
      y(step+1) = y(step) + sign(randn);
      % plot(x, y, 'b');
      % grid on
      % hold on
    end
    R2(i, samples) =  (x(end)^2+y(end)^2)^2;
  end
  % hold off
end
subplot(2,1,1)
plot(sqrt(mean(R2,2)), sqrt(N_incr))
title('Distance vs sample size')
grid on
xlabel('N^{1/2}')
ylabel('<R^2>^{1/2}')
subplot(2,1,2)
loglog(sqrt(mean(R2,2)), sqrt(N_incr))
title('Loglog distance vs sample size')
xlabel('N^{1/2}')
ylabel('<R^2>^{1/2}')
grid on
