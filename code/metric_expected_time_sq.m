function distance = metric_expected_time_sq(dt, cry_vec, cry_vec0, arrival_time, integral_start, integral_end)

% by default, integration starts from arrival time, ends at 600

% we calculate \frac{\int_{arrival_time}^{end_time} t [f(t) - f_0(t)]^2 dt}{\int_{arrival_time}^{end_time} [f(t) - f_0(t)]^2 dt}


% f(t) is the cry_vec

% f_0(t) is the cry_vec0, where you are originally at the destination time
% zone


denom = metric_total_variation_square(dt, cry_vec, cry_vec0, integral_start, integral_end);

sq_variation_vec = (cry_vec - cry_vec0).^2;
numer = 0;
for k = integral_start:integral_end
    numer = numer + dt * sq_variation_vec(k)*(dt*k-arrival_time);
end

distance = numer / denom;
end