% defining the given analytic sltn
function y = analytic(a, r)

    y = ((1/(1-a^2))^2)*(exp(-a*r) - exp(-r).*(1+0.5*((1-a^2).*r))); 

end 
