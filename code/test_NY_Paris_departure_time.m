close all
clear all
clc

dd = 0.005;
tvec = 0:dd:800;


% experiment
morning = 242; % 8am
noon = 240 + 6; % 12 at noon
evening = 240 + 6 + 8; % 8 pm

flight_time = 7;
time_difference = 24-6; 

plot_start = round(151/dd);
plot_end = round(500/dd);

integral_start = round((noon+flight_time)/dd)+1;
integral_end = round(600/dd);

%{
INPUTS:
with_light:         1 if light is 12-12 alternating based on your time zone,
                    0 if constant 24hour light

with_food:          1 if food is regularly served at 8, 12, and 19 everyday, based on your time zone
                    0 if food is held constant 

departure_time:     the time your flight begins, PLEASE choose a number > 120 to allow the system settles down first

arrival_time:       the time your flight arrives, = departure_time + flight_time

time_difference:    the number of time zones between your origin and your destination, counting from west to east

meal_time_on_plane: the time of meals you want to specify on plane, an array whose ...
                    length is equal to the number of meals on the plane. every number ...
                    in the array must be betweeen departure_time and arrival_time

meal_ampl_on_plane: of the same length as meal_time_on_plane


%}
% always stay at origin
[ampk_vec_ori, lighting_per_max_vec_ori, lighting_cry_max_vec_ori,...
    lighting_rev_max_vec_ori, lighting_ror_max_vec_ori, per_vec_ori, ...
    cry_vec_ori, rev_vec_ori, ror_vec_ori, bmal_vec_ori] = eating_circadian(dd,1,1,800,800,0,[],[]);

% always stay at destination
[ampk_vec_dest, lighting_per_max_vec_dest, lighting_cry_max_vec_dest,...
    lighting_rev_max_vec_dest, lighting_ror_max_vec_dest, per_vec_dest, ...
    cry_vec_dest, rev_vec_dest, ror_vec_dest, bmal_vec_dest] = eating_circadian(dd,1,1, 0,0,time_difference,[],[]);


[ampk_vec_1, lighting_per_max_vec_1, lighting_cry_max_vec_1,...
    lighting_rev_max_vec_1, lighting_ror_max_vec_1, per_vec_1, ...
    cry_vec_1, rev_vec_1, ror_vec_1, bmal_vec_1] = eating_circadian(dd,1,1, evening,evening + flight_time,time_difference,[],[]);

[metric_total_variation_square(dd, cry_vec_1, cry_vec_dest, integral_start, integral_end), ...
    metric_total_variation_abs(dd, cry_vec_1, cry_vec_dest, integral_start, integral_end),...
    metric_expected_time_sq(dd, cry_vec_1, cry_vec_dest, evening+flight_time, integral_start, integral_end),...
    metric_expected_time_abs(dd, cry_vec_1, cry_vec_dest, evening+flight_time, integral_start, integral_end)]
%{
ampk_mat = zeros(length(ampk_vec_ori), 15);
lighting_per_max_mat = zeros(length(ampk_vec_ori),15);
lighting_cry_max_mat = zeros(length(ampk_vec_ori),15);
lighting_rev_max_mat = zeros(length(ampk_vec_ori),15);
lighting_ror_max_mat = zeros(length(ampk_vec_ori),15);
per_mat = zeros(length(ampk_vec_ori),15);
cry_mat = zeros(length(ampk_vec_ori),15);
rev_mat = zeros(length(ampk_vec_ori),15);
ror_mat = zeros(length(ampk_vec_ori),15);
bmal_mat = zeros(length(ampk_vec_ori),15);

sq_metric_vec = zeros(15,1);
abs_metric_vec = zeros(15,1);
Et_metric_sq_vec = zeros(15,1);
Et_metric_abs_vec = zeros(15,1);


% see when shall we give a meal
for count = 2:15
    [ampk_mat(:,count),lighting_per_max_mat(:,count),lighting_cry_max_mat(:,count),...
        lighting_rev_max_mat(:,count),lighting_ror_max_mat(:,count), per_mat(:,count),...
        cry_mat(:,count), rev_mat(:,count), ror_mat(:,count), bmal_mat(:,count)]...
        = eating_circadian(dd,1,1,noon,noon+flight_time,time_difference,[noon+count-1],[1]);
    
    sq_metric_vec(count) = metric_total_variation_square(dd, cry_mat(:,count), cry_vec_dest, integral_start, integral_end);
    abs_metric_vec(count)= metric_total_variation_abs(dd, cry_mat(:,count), cry_vec_dest, integral_start, integral_end);
    Et_metric_sq_vec(count) = metric_expected_time_sq(dd, cry_mat(:,count), cry_vec_dest, noon+flight_time, integral_start, integral_end);
    Et_metric_abs_vec(count) = metric_expected_time_abs(dd, cry_mat(:,count), cry_vec_dest, noon+flight_time, integral_start, integral_end);
    
end
[ampk_mat(:,1),lighting_per_max_mat(:,1),lighting_cry_max_mat(:,1),...
        lighting_rev_max_mat(:,1),lighting_ror_max_mat(:,1), per_mat(:,1),...
        cry_mat(:,1), rev_mat(:,1), ror_mat(:,1), bmal_mat(:,1)]...
        = eating_circadian(dd,1,1,noon,noon+flight_time,time_difference,[],[]);
sq_metric_vec(1) = metric_total_variation_square(dd, cry_mat(:,1), cry_vec_dest, integral_start, integral_end);
abs_metric_vec(1)= metric_total_variation_abs(dd, cry_mat(:,1), cry_vec_dest, integral_start, integral_end);
Et_metric_sq_vec(1) = metric_expected_time_sq(dd, cry_mat(:,1), cry_vec_dest, noon+flight_time, integral_start, integral_end);
Et_metric_abs_vec(1) = metric_expected_time_abs(dd, cry_mat(:,1), cry_vec_dest, noon+flight_time, integral_start, integral_end);
hold on
plot([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14],sq_metric_vec, '-*','LineWidth',2,'MarkerSize',6)
plot([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14],abs_metric_vec,'-*','LineWidth',2,'MarkerSize',6)
plot([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14],Et_metric_sq_vec, '-*','LineWidth',2,'MarkerSize',6)
plot([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14],Et_metric_abs_vec, '-*','LineWidth',2,'MarkerSize',6)
legend('sq','abs','Et, weight based on sq', 'Et, weight based on abs')

figure(2)
%subplot(2,1,1)
%hold on
%plot(tvec(plot_start:plot_end), lighting_rev_max_mat(plot_start:plot_end,1), '-*', 'MarkerSize', 1)
%plot(tvec(plot_start:plot_end), lighting_rev_max_vec_ori(plot_start:plot_end))
%plot(tvec(plot_start:plot_end), lighting_rev_max_vec_dest(plot_start:plot_end))
%plot(tvec(plot_start:plot_end), 1-ampk_mat(plot_start:plot_end,1)) % no meal
%plot(tvec(plot_start:plot_end), 1-ampk_mat(plot_start:plot_end,2)) % meal at t=1
%plot(tvec(plot_start:plot_end), 1-ampk_mat(plot_start:plot_end,6)) % meal at t=5
%legend('time on your phone', 'time at origin' ,'time at destination')

%subplot(2,1,2)
hold on
plot(tvec(plot_start:plot_end), lighting_rev_max_mat(plot_start:plot_end,1))

plot(tvec(plot_start:plot_end), cry_mat(plot_start:plot_end,1), '-*', 'MarkerSize', 1)
plot(tvec(plot_start:plot_end), cry_mat(plot_start:plot_end,2), '-*', 'MarkerSize', 1)
plot(tvec(plot_start:plot_end), cry_mat(plot_start:plot_end,6), '-*', 'MarkerSize', 1.4)

plot(tvec(plot_start:plot_end), cry_vec_ori(plot_start:plot_end), 'LineWidth', 0.8)
plot(tvec(plot_start:plot_end), cry_vec_dest(plot_start:plot_end), 'LineWidth', 0.8)

legend('time on your phone', 'no meal','meal at t=1', 'meal at t=5', 'triplet at origin', 'triplet at destination')


set(gcf,'Position',[200 200 1000 600])
%}



