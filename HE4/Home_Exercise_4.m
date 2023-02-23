clear all 
% Samples and walks
N = 50000;
N_tests = 10000; % number of times to test each N-value 
N_incr = linspace(1000,N,N/10000); % test N-values 

test_results = zeros(1, length(N_incr)); 



for i = 1:length(N_incr)
   N_test = N_incr(i); 
   
    for j = 1:N_tests
    % each test starts from the origin 
      position = [0,0]; 
        for k = 1:N_test
            
            position = position + random_step(); 

        end
        
      test_results(i) = test_results(i) + sqrt(position(1)^2 + position(2)^2);
      
    end 

end

test_results = test_results./N_tests; % averages the distances 

hold on 
plot(test_results, sqrt(N_incr), 'rs-')
title('Distance vs sample size')
grid on
xlabel('N^{1/2}')
ylabel('<R^2>^{1/2}')

hold on

% function generates a random vector
% is either [0,0], [1,0], [-1,0], [0, 1], [0, -1]
function step = random_step()

possibles = [[0,0]; [1,0]; [-1,0]; [0, 1]; [0, -1]]; 

step = possibles(randi(5), :); % function returns a random step 

end 