% defining a numeric solution to the entire problem 

function y = numeric_phi(r, a, r_max)

N  = 1000; % number of integration steps 

% we split the sltn in two for readability 
first = phi_great(a, r)*bode(@less_integ, 0, r, N, a);  
second = phi_less(a, r)*bode(@great_integ, r, r_max, N, a); 

y = first + second;  

end 
