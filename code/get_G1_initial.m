function result = get_G1_initial(t0)


% default tau = 1.4
tau = 1.4;

% we have G1(t0) = tau^2 * f''(t0) + 2*tau*f'(t0) + f(t0)
% second order centered finite difference for f''
f_prime = (meals(t0 + 0.0001)-meals(t0-0.0001))/(0.0002);
f_primeprime = (meals(t0+0.0001)-2*meals(t0)+meals(t0-0.0001))/(0.0001^2);
result = tau^2 * f_primeprime + 2*tau*f_prime + meals(t0);
end