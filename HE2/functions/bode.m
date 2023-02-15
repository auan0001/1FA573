% numerical integration using Bode's method 
% Tests can be found in lecture_1_2.m

function integ = bode(f, a, b, N, a_val)

  % f - function handle, to be integrated 
  % a , b - upper and lower bounds 
  % N, number of steps 
  % a_val - given a-constant 

  % assuming N is a multiple of 4

  h = (b-a)/N; % step-length

  integral = 0; % prealocation 

  for n = 1:4:N-3

    part_1 = 7*f(a_val, a+(n-1)*h) + 32*f(a_val, a+(h*n)); 
    part_2 = 12*f(a_val, a+h*(n+1)) + 32*f(a_val, a+h*(n+2)) + 7*f(a_val, a+h*(n+3)); 
    integral = integral + (((2*h)/45) * (part_1 + part_2)); 

  end

  integ = integral; %returning the final value 

end 
