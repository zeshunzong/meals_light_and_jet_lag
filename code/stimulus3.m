function stimulus_vec = stimulus3(timevec, dt, meal_time, meal_amplitude)
% amplitude_vec is a vector of length(timevec)
% for each clock = 1:length(timevec), stimulus_vec(clock) = 4.1 *
% amplitude_of_that_meal * delta(clock*dt- mealTime)
stimulus_vec = zeros(length(timevec), 1);
pointer = 1;
for clock = 1:length(timevec)
    % currentTime = clock*dt
    % if at currentTime we have a meal, there should be a corresponding
    % amplitude in stimulus_vec(clock)
    currentTime = clock*dt;
    % check if currentTime is within the next meal block
    if pointer <= length(meal_time) && currentTime >= meal_time(pointer)-0.5*dt && currentTime < meal_time(pointer)+0.5*dt
        stimulus_vec(clock) = 4.1* meal_amplitude(pointer) * 1/dt;
        pointer = pointer + 1;
    end
end
end