function epsilon_vec = lighting_through_epsilon(departure_time, arrival_time, endtime, dt, timevec, epsilon0, time_difference)

% time difference is the number of time zones between the two district,
% counting from west to east. In this way, the destination is
% time_difference hours earlier (in the sense that the time has passed) the
% origin.

% daylight at origin: ZT0-ZT12
% daylight at destination: ZT(x)-ZT(12+x mod 24)

% before departure, we give regular alternative lighting, starting with 12h light,
% then 12h dark, then 12h light, etc, assuming we are in Shanghai time zone

epsilon_vec = zeros(length(timevec)+round(24/dt), 1);
for clock = 1 : length(timevec)+round(24/dt)
    remainder = mod(clock*dt, 24);
    if remainder <= 12
        epsilon_vec(clock) = epsilon0;
    else
        epsilon_vec(clock) = -epsilon0;
    end
end

% during the flight, can be entered manually
for clock = round(departure_time/dt) + 1 : round(arrival_time / dt)
    % assuming darkness
    epsilon_vec(clock) = -epsilon0;
end

% after arriving at destination 
%morning_start = 0 + time_difference;
%evening_start = 12 + time_difference; % this may be larger than 24

for clock = round(arrival_time / dt) + 1 : round(endtime/dt)
    
    % the sunshine at destination would be the same as the sunshine
    % time_difference hours later at origin
    epsilon_vec(clock) = epsilon_vec(clock + round(time_difference / dt));
end

epsilon_vec = epsilon_vec(1:length(timevec));

end



    
   