% RK2-function which returns the function value at the next time-incident 
function res = RK2(func, y_calc, y_cur, t, h)

    k = h*func(t, y_calc); 
    res = y_cur + h*func(t + 0.5*h, y_calc + 0.5*k); 

end 
