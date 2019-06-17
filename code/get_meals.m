function [meal_time, meal_amplitude] = get_meals(departure_time, arrival_time, endtime, time_difference, meal_time_on_plane, meal_ampl_on_plane)

% return two vectors, indicating time of meals and amplitude of meals,
% before departure_time the meals are specified at ZT2, ZT6, and ZT13, with amplitudes
% equal to 0.6, 1, 1 respectively

meal_time = [2, 6, 13];
meal_amplitude = [0.6, 1, 1];
i=4;
while 1
    temp = meal_time(i-3)+24;
    if temp > departure_time
        break
    end
    meal_time(i) = temp;
    meal_amplitude(i) = meal_amplitude(i-3);
    i = i + 1;
end
% 252 is good

% for meals between departure time and arrival time, i.e. during the flight,
% they are specified manually. 
meal_time_during_flight = meal_time_on_plane;
meal_amplitude_during_flight = meal_ampl_on_plane;
% notice that every number in the vector meal_time_during_flight should be
% in (departure_time, arrival_time).
meal_time = [meal_time, meal_time_during_flight];
meal_amplitude = [meal_amplitude, meal_amplitude_during_flight];

%{
if food_adjust == 5

    % for meals between departure time and arrival time, i.e. during the flight,
    % they are specified manually. 

    meal_time_during_flight = [258];
    meal_amplitude_during_flight = [1];
    % notice that every number in the vector meal_time_during_flight should be
    % in (departure_time, arrival_time).
    meal_time = [meal_time, meal_time_during_flight];
    meal_amplitude = [meal_amplitude, meal_amplitude_during_flight];
    %meal_time
end
%}

%%%%%%%%%%%
% after flight

if arrival_time == 0
    i = 1;
    currentTime=0;
else
    i = length(meal_time)+1;
    currentTime = max([arrival_time, 14]);
end

breakfast_time = mod(2 + (24-time_difference), 24);
lunch_time = mod(6 + (24-time_difference), 24);
dinner_time = mod(13 + (24-time_difference), 24);

while currentTime <= endtime
    if mod(currentTime, 24) == breakfast_time
        meal_time(i) = currentTime;
        meal_amplitude(i) = 0.6;
        i = i + 1;
    elseif mod(currentTime, 24) == lunch_time
        meal_time(i) = currentTime;
        meal_amplitude(i) = 1;
        i = i + 1;
    elseif mod(currentTime, 24) == dinner_time
        meal_time(i) = currentTime;
        meal_amplitude(i) = 1;
        i = i + 1;
    end
    currentTime = currentTime + 1;
end
if meal_time(1) == 0
    meal_time = meal_time(2:end);
    meal_amplitude = meal_amplitude(2:end);
end

