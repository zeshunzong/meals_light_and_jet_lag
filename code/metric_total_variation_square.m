function distance = metric_total_variation_square(dt, cry_vec, cry_vec0, integral_start, integral_end)

% by default, integration starts from arrival time, ends at 600

% we calculate \int_{arrival_time}^{end_time} [f(t) - f_0(t)]^2 dt

% f(t) is the cry_vec
% f_0(t) is the cry_vec0, where you are originally at the destination time
% zone

sq_variation_vec = (cry_vec - cry_vec0).^2;


distance = 0;
for k = integral_start:integral_end
    distance = distance + dt * sq_variation_vec(k);
end

end