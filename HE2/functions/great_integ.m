% defining integrand for phi_great term 
function y = great_integ(a, r)

source = -((r*exp(-r))/2); 
y = phi_great(a, r) * source; 

end 
