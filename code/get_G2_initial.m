function result = get_G2_initial(t0)

% return the initial value of G2 at t0, where t0 is the starting time 
% from which we deviate from the regular schedule

% simply use finite difference for f'
% recall that G2(t0) = tau * f'(t0) + f(t0)

tau = 1.4;

f_prime = (meals(t0 + 0.0001)-meals(t0-0.0001))/(0.0002);
result = tau*f_prime + meals(t0);


end