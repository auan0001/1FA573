% RK3-function which returns the function value at the next time-incident 
function res = RK3(func, y_calc, y_cur, t, h)

    k1 = h*func(t, y_calc); 
    k2 = h*func(t+h*(0.5), y_calc + k1/2);
    k3 = h*func(t+h, y_calc - k1 + 2*k2);
    res = y_cur + (1/6)*(k1 + 4*k2 + k3);

end 
