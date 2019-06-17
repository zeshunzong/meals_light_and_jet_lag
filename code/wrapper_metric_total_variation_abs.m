function distance = wrapper_metric_total_variation_abs(dt, with_light, with_food,...
    departure_time, arrival_time, time_difference, meal_time_on_plane,...
    meal_ampl_on_plane)

% by default, integration starts from arrival time, ends at 600

% we calculate \int_{arrival_time}^{end_time} |f(t) - f_0(t)| dt


% f(t) is the cry_vec
[ampk_vec, lighting_per_max_vec, lighting_cry_max_vec,...
    lighting_rev_max_vec, lighting_ror_max_vec, per_vec, ...
    cry_vec, rev_vec, ror_vec, bmal_vec] = eating_circadian(dt, with_light, with_food,...
    departure_time, arrival_time, time_difference, meal_time_on_plane,...
    meal_ampl_on_plane);

% f_0(t) is the cry_vec0, where you are originally at the destination time
% zone
[ampk_vec0, lighting_per_max_vec0, lighting_cry_max_vec0,...
    lighting_rev_max_vec0, lighting_ror_max_vec0, per_vec0, ...
    cry_vec0, rev_vec0, ror_vec0, bmal_vec0] = eating_circadian(dt, with_light, with_food,...
    0, 0, time_difference, meal_time_on_plane,...
    meal_ampl_on_plane);

abs_variation_vec = abs(cry_vec - cry_vec0);

integral_start = round(arrival_time/dt) + 1;
integral_end = round(600/dt);

distance = 0;
for k = integral_start:integral_end
    distance = distance + dt * abs_variation_vec(k);
end

end