function ampk_vec = AMPK_of_t_profile(departure_time, arrival_time, endtime, dt, timevec, time_difference, meal_time_on_plane, meal_ampl_on_plane)

% t0 is the starting time 
% from which we deviate from the regular schedule
% before t0, meals are supplied regularly
% after t0, meals are determined by our stimulus function
% we are solving ODEs numerically, so only at discrete 
% t = clock*dt

[meal_time, meal_amplitude] = get_meals(departure_time, arrival_time, endtime, time_difference, meal_time_on_plane, meal_ampl_on_plane);

% values of the ampk_vec should be determined by
% solving the ODEs for G1, G2 and G3




sti_vec = stimulus3(timevec, dt, meal_time, meal_amplitude);

G1_vec = zeros(length(timevec), 1);
G2_vec = zeros(length(timevec), 1);
G3_vec = zeros(length(timevec), 1);
G1 = get_G1_initial(0);
G2 = get_G2_initial(0);
G3 = meals(0);

tau = 1.4;
ampk_vec = zeros(length(timevec), 1);

for clock = 1: length(timevec)
    G1_new = (sti_vec(clock)-G1)*dt/tau+G1;
    
    G2_new = (G1_new - G2)*dt/tau + G2;
    
    G3_new = (G2_new - G3)*dt/tau + G3;
    
    G1_vec(clock) = G1_new;
    G2_vec(clock) = G2_new;
    G3_vec(clock) = G3_new;
    G1 = G1_new;
    G2 = G2_new;
    G3 = G3_new;
    
    % since G3 = f
    % we also want it to be bounded between 0 and 1
    if 1-G3_new < 0.000001
        ampk_vec(clock) = 0.00001;
        G3 = 0.99999;
    elseif 1-G3_new > 0.99999
        ampk_vec(clock) = 0.99999;
        G3 = 0.00001;
    else
        ampk_vec(clock) = 1 - G3_new;
    end
end
end




     