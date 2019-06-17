function distance = metric_total_variation_abs(dt, cry_vec, cry_vec0, integral_start, integral_end)

% by default, integration starts from arrival time, ends at 600

% we calculate \int_{arrival_time}^{end_time} |f(t) - f_0(t)| dt


% f(t) is the cry_vec
% f_0(t) is the cry_vec0, where you are originally at the destination time
% zone

abs_variation_vec = abs(cry_vec - cry_vec0);
distance = 0;

for k = integral_start:integral_end
    distance = distance + dt * abs_variation_vec(k);
end

end