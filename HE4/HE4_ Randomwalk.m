clear all 
% Samples and walks
N = 5000;
N_tests = 1000; % number of times to test each N-value 
N_incr = linspace(1000,N,N/10000); % test N-values 

test_results = zeros(1, length(N_incr)); 

% initialize the position vector
pos = [0,0];

% create an array to store the position at each step
pos_array = zeros(N,2);

for k = 1:N
    % generate a random step and update the position
    step = random_step();
    pos = pos + step;
    
    % store the position in the array
    pos_array(k,:) = pos;
end

% plot the travel
plot(pos_array(:,1),pos_array(:,2),'b-')
title('Example of a random walk')
xlabel('x-axis')
ylabel('y-axis')

% function generates a random vector
% is either [0,0], [1,0], [-1,0], [0, 1], [0, -1]
function step = random_step()

possibles = [[0,0]; [1,0]; [-1,0]; [0, 1]; [0, -1]]; 

step = possibles(randi(5), :); % function returns a random step 

end 
