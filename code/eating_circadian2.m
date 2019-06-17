function [ampk_vec, lighting_per_max_vec, lighting_cry_max_vec,...
    lighting_rev_max_vec, lighting_ror_max_vec, per_vec, ...
    cry_vec, rev_vec, ror_vec, bmal_vec] = eating_circadian2(...
    dt, with_light, with_food, departure_time, arrival_time, time_difference, meal_time_on_plane, meal_ampl_on_plane)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
INPUTS:
dt:                 timestep
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
The following are the variables:
per         Concentration of Per mRNA
cry         Concentration of Cry mRNA
rev         Concentration of Rev-Erb mRNA
ror         Concentration of Ror mRNA
bmal        Concentration of Bmall mRNA
PER         Concentration of PER protein
CRY         Concentration of CRY protein
REV         Concentration of REV-ERB protein
ROR         Concentration of ROR protein
BMAL        Concentration of BMAL1 protein
PC          Concentration of PER-CRY 
CB          Concentration of CLOCK-BMAL1 protein complex
nampt       Concentration of Nampt mRNA
NAMPT       Concentration of NAMPT protein
dbp         Concentratino of Dbp mRNA
NAD         NAD level
Act_SIRT    activity of SIRT1
Act_PGC1a   activity of PGC1a
Atc_AMPK    activity of AMPK
PGC1a       nuclear abundance of PGC1a
agonist_rev Concentration of REV-ERB agonist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Two expressions:
Ki_bmal_rev threshold of repression of Bmal by Rev-Erb
Ki_cry_rev  threshold of repression of cry by Rev-Erb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The following are the parameters:

Degradation rates, unit 1/h (15 parameters)
dm_bmal         Bmal mRNA degradation rate
dm_cry          Cry mRNA degradation rate
dm_dbp          Dbp mRNA degradation rate
dm_nampt        Nampt mRNA degradation rate
dm_per          Per mRNA degradation rate
dm_rev          Rev-Erb mRNA degradation rate
dm_ror          Ror mRNA degradation rate
dp_cry          CRY protein degradation rate
dp_nampt        NAMPT protein degradation rate
dp_per          PER protein degradation rate
dp_rev          REV-ERB protein degradation rate
dp_ror          ROR degradation rate
d_cb            CLOCK-BMAL complex degradation rate
d_pc            PER-CRY complex degradation rate

Complexation kinetics rates (4 parameters)
kass_cb         CLOCK-BMAL association rate l/(nmol*h)
kass_pc         PER-CRY association rate l/(nmol*h)
kdiss_cb        CLOCK-BMAL dissociation rate 1/h
kdiss_pc        PER-CRY dissociation rate 1/h

Maximal transcription rates, unit nmol/(l*h) (7 parameters)
Vmax_bmal       Bmal maximal transcription rate
Vmax_cry        Cry maximal transcription rate
Vmax_dbp        Dbp maximal transcription rate
Vmax_nampt      Nampt maximal transcription rate
Vmax_per        Per maximal transcription rate
Vmax_rev        Rev-Erb maximal transcription rate
Vmax_ror        Ror maximal transcription rate

Activation ratios, dimensionless (7 parameters)
fold_bmal       Activation ratio of Bmal by ROR
fold_cry        Activation ratio of Cry by CLOCK-BMAL
fold_dbp        Activation ratio of Dbp by CLOCK-BMAL
fold_nampt      Activation ratio of Nampt by CLOCK-BMAL
fold_per        Activation ratio of Per by CLOCK-BMAL
fold_rev        Activation ratio of Rev-Erb by CLOCK-BMAL
fold_ror        Activation ratio of Ror by CLOCK-BMAL

Regulation thresholds, unit l/nmol (15 parameters)
Ka_bmal_ror     Regulation threshold of Bmal by ROR
Ka_cry_cb       Regulation threshold of Cry by CLOCL-BMAL
Ka_dbp_cb       Regulation threshold of Dbp by CLOCL-BMAL
Ka_nampt_cb     Regulation threshold of Nampt by CLOCL-BMAL
Ka_per_cb       Regulation threshold of Per by CLOCL-BMAL
Ka_rev_cb       Regulation threshold of Rev-Erb by CLOCL-BMAL
Ka_ror_cb       Regulation threshold of Ror by CLOCL-BMAL
Ki_bmal_rev0    Regulation threshold of Bmal by Rev-Erb
Ki_cry_rev0     Regulation threshold of Cry by Rev-Erb
Ki_cry_pc       Regulation threshold of Cry by PER-CRY
Ki_dbp_pc       Regulation threshold of Dbp by PER-CRY
Ki_nampt_pc     Regulation threshold of Nampt by PER-CRY
Ki_per_pc       Regulation threshold of Per by PER-CRY
Ki_rev_pc       Regulation threshold of Rev-Erb by PER-CRY
Ki_ror_pc       Regulation threshold of Ror by PER-CRY

Hill coefficients

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% asssign parameter values
% degradation rates
[dm_bmal,dm_cry,dm_dbp,dm_nampt,dm_per,dm_rev,dm_ror,dp_bmal, ...
    dp_cry,dp_nampt,dp_per,dp_rev,dp_ror,d_cb,d_pc] = deal(0.827333442085,...
    0.319848706181,0.384486532062, ...
    0.814311309051,0.323114598647,4.0479072846,0.246575760727,...
    0.186413466797,0.599026119971,49.8841023982,10.9446515392,...
    0.281564063057,0.0340112196281,0.197714012552,0.609290138329);

% complexation kinetic rates
[kass_cb, kass_pc, kdiss_cb, kdiss_pc] = deal(0.0162923358655,12.302005485,...
    0.00613502224231, 0.0365867175408);

% maximal transcription rates
[Vmax_bmal, Vmax_cry, Vmax_dbp, Vmax_nampt, Vmax_per, Vmax_rev, Vmax_ror]...
    = deal(0.0916862770193, 0.702216664807, 0.0802009609453,...
    3.49035201578, 0.730201742662, 1.12297601784, 6.9843472736);

% activation ratios
[fold_bmal,fold_cry,fold_dbp,fold_nampt,fold_per,fold_rev,fold_ror] = ...
    deal(15.9258093373, 1.1604489571, 400.0, 1.57880681573, ...
    12.977351123, 73.2342431701, 335.923333883);

% regulation thresholds
[Ka_bmal_ror, Ka_cry_cb, Ka_dbp_cb, Ka_nampt_cb, Ka_per_cb, Ka_rev_cb,...
    Ka_ror_cb, Ki_bmal_rev0, Ki_cry_rev0, Ki_cry_pc, Ki_dbp_pc, ...
    Ki_nampt_pc, Ki_per_pc, Ki_rev_pc, Ki_ror_pc] = ...
    deal(0.00560038941378, 1.0089387144, 0.308745016237, 3.54586790835,...
    2.03485134763, 0.260846828116, 0.266407416327, 0.0108449480001, ...
    0.248955507809, 0.00338463577329, 2.23913672671, 0.0137106537972, ...
    0.273493946059, 28.5885406354, 0.0072858432208);

% Hill coefficients
[hill_bmal_rev, hill_bmal_ror, hill_cry_cb, hill_cry_pc, hill_cry_rev, ...
    hill_dbp_cb, hill_dbp_pc, hill_nampt_cb, hill_nampt_pc, hill_per_cb,...
    hill_per_pc, hill_rev_cb, hill_rev_pc, hill_ror_cb, hill_ror_pc] = ...
    deal(4.32985205032, 1.83992599578, 9.1109447538, 2.43715119318, ...
    4.20952050286, 7.32066818222, 10.4312927466, 1.91470474775, ...
    1.34080593157, 8.52414053707, 8.53897990872, 9.83701536835, ...
    3.31257899336, 9.36456505724, 1.84102439743);

% translation rates
[kp_bmal, kp_cry, kp_nampt, kp_per, kp_rev, kp_ror] = ...
    deal(0.628507384997, 3.76912711677, 58.9647983852, 13.2909782781, ...
    0.0513221194724, 0.0412765888526);

% protein stability modulation constants
[m_cry_ampk, m_nampt_ampk, m_per_ampk, m_per_sirt] = ...
    deal(0.07940386211, 0.6352243791, 0.005243953731, 0.005452322816);

% NAD kinetics, Sirt1 and PGC1a activity
[Vsirt,Ksirt,d_nad,Knad,NAD_basal,Vnad,NAD_tot,Knam,Vpg,Kpg1,Kpg2] = ...
    deal(0.915854846427, 0.75, 378.009130749, 0.321746039086, ...
    0.9116166306, 296.3933869, 4.166901679, 2.76496, 24.06372088, ...
    0.046630145542, 12.3526351747);

% Pulse parameters
[tc1, tc2, tc3, Td1, Td2, Td3, Tr1, Tr2, Tr3, amp1, amp2, amp3] = ...
    deal(4.38900149, 15.75, 18.875, 2.25, 1.5, 15.25, 2.6, 1.8, ...
    0.5, 6.0, 0.9778008, 0.803062);

% Chronotherapy timings
[tc4, Td4, Tr4, amp4] = deal(13.664, 2.83718, 1.86794, 0.465852);

Csirt = 1;
Campk = 1;
offs = 0.02;
%{
% Miscellaneous constants used to describe perturbations
status = 'WT';

switch status
    case 'WT'
        Csirt = 1;
        Campk = 1;
        offs = 0.02;
    case 'SIRT1 KO'
        Csirt = 0;
        Campk = 1;
        offs = 0.02;
    case 'LKB1 KO'
        Csirt = 1;
        Campk = 0.0375;
        offs = 0.02;
    case 'HFD'
        Csirt = 1;
        Campk = 0.05;
        offs = 0.02;
    case 'fasting'
        Csirt = 1;
        Campk = 0.05;
        offs = 2.6;
end
%}


[per, cry, rev, ror, bmal, PER, CRY, REV, ROR, BMAL, PC, CB, nampt, NAMPT,...
    dbp, NAD, Act_SIRT, Act_PGC1a, Act_AMPK, PGC1a, agonist_rev, ...
    Ki_bmal_rev, Ki_cry_rev] = ...
    deal(0.163324778297151, 1.575830320129322, 0.282259161402967, 0.227099482727232,...
    0.766426221859058, 0.198679063661776, 1.950748229899452, 0.053675150905969,...
    0.226476576094958, 1.947963174964374, 6.959449010157386, 0.168452286394014,...
    0.153705591476109, 0.222553862031172, 0.007927828871210, 0.945227127405072,...
    0.510738763256038, 0.290812134965029, 0.178438030440945, 0.803061999996042,...
    0.002302193822186, 0.010820011696557, 0.248383072596635);



per = 0.107693115132316;
cry = 2.018311284858204;
rev = 0.195991369800841;
ror = 0.837184709330401;
bmal = 1;
PER = 1;
CRY = 1;
REV = 1;
ROR = 1;
BMAL = .2;
PC = .2;
CB = .2;
nampt = .2;
NAMPT = .2;
dbp = .2;
NAD = .2;
Act_SIRT = .2;
Act_PGC1a = .2;
Act_AMPK = .2;
PGC1a = .2;
agonist_rev = .2;
Ki_bmal_rev = .2;
Ki_cry_rev = .2;



% expression of changing max transcription rate of per
lighting_per_max = Vmax_per; % initialized with no lighting effect
lighting_cry_max = Vmax_cry;
lighting_ror_max = Vmax_ror;
lighting_rev_max = Vmax_rev;

% initialize time variables
% unit: hours

endtime = 800;
timevec = 0:dt:endtime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT into the model
%departure_time = 246;
%departure_time = 6+24*5;
%arrival_time = 246+15;
%arrival_time = 6+24*5+5;
%time_difference = 12;
%time_difference = 21;
% food intake
ampk_vec = AMPK_of_t_profile(departure_time, arrival_time, endtime, dt, timevec, time_difference, meal_time_on_plane, meal_ampl_on_plane);
% light
epsilon0 = 0.3;
epsilon_vec = lighting_through_epsilon(departure_time, arrival_time, endtime, dt, timevec, epsilon0, time_difference);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% let the plot start from 3/4 of the timevec, avoid initial states
plot_start = floor(1/dt);
plot_end = floor(600/dt);

% vectors to store values at each iteration 
per_vec = zeros(length(timevec), 1);
cry_vec = zeros(length(timevec), 1);
rev_vec = zeros(length(timevec), 1);
ror_vec = zeros(length(timevec), 1);
bmal_vec = zeros(length(timevec), 1);
PER_vec = zeros(length(timevec), 1);
CRY_vec = zeros(length(timevec), 1);
REV_vec = zeros(length(timevec), 1);
ROR_vec = zeros(length(timevec), 1);
BMAL_vec = zeros(length(timevec), 1);
PC_vec = zeros(length(timevec), 1);
CB_vec = zeros(length(timevec), 1);
nampt_vec = zeros(length(timevec), 1);
NAMPT_vec = zeros(length(timevec), 1);
dbp_vec = zeros(length(timevec), 1);
NAD_vec = zeros(length(timevec), 1);
Act_SIRT_vec = zeros(length(timevec), 1);
Act_PGC1a_vec = zeros(length(timevec), 1);
Act_AMPK_vec = zeros(length(timevec), 1);
PGC1a_vec = zeros(length(timevec), 1);
agonist_rev_vec = zeros(length(timevec), 1);
% expressions
Ki_bmal_rev_vec = zeros(length(timevec), 1);
Ki_cry_rev_vec = zeros(length(timevec), 1);

% additionally, Vmax_per can be changed
lighting_per_max_vec = zeros(length(timevec), 1);
lighting_cry_max_vec = zeros(length(timevec), 1);
lighting_ror_max_vec = zeros(length(timevec), 1);
lighting_rev_max_vec = zeros(length(timevec), 1);


% a switch on food
if with_food == 1
    
elseif with_food == 0
    ampk_vec = ones(length(ampk_vec),1).*0.5;
end

% time loop

for clock = 1:length(timevec)
    
    % a total of 23 equations
    
    % 1
    numer = lighting_per_max * (1+fold_per*(CB/(Ka_per_cb*(1+Act_SIRT)))^hill_per_cb);
    denom = 1 + (CB/Ka_per_cb/(1+Act_SIRT))^hill_per_cb * (1 + (PC/Ki_per_pc)^hill_per_pc);
    per_new = per + ((-dm_per * per) + numer/denom) * dt;
  
    % 2
    numer = lighting_cry_max *(1+fold_cry*(CB/Ka_cry_cb/(1+Act_SIRT))^hill_cry_cb);
    denom = (1+(CB/Ka_cry_cb/(1+Act_SIRT))^hill_cry_cb * (1 + (PC/Ki_cry_pc)^hill_cry_pc))*(1+(REV/Ki_cry_rev)^hill_cry_rev);
    cry_new = cry + (-dm_cry * cry + numer/denom) * dt;

    % 3
    numer = lighting_rev_max *(1 + fold_rev *(CB/Ka_rev_cb/(1+Act_SIRT))^hill_rev_cb);
    denom = 1 + (CB/Ka_rev_cb/(1+Act_SIRT))^hill_cry_cb * (1+ (PC/Ki_rev_pc)^hill_rev_pc);
    rev_new = rev + (-dm_rev*rev + numer/denom) * dt;
    
    % 4
    numer = lighting_ror_max * (1 + fold_ror * (CB/Ka_ror_cb/(1+Act_SIRT))^hill_ror_cb);
    denom = 1 + (CB/Ka_ror_cb/(1+Act_SIRT))^hill_ror_cb * (1 + (PC/Ki_ror_pc)^hill_ror_pc);
    ror_new = ror + (-dm_ror * ror + numer/denom) * dt;

    % 5
    numer = Vmax_bmal * (1 + fold_bmal * (1+Act_PGC1a) * (ROR/Ka_bmal_ror)^hill_bmal_ror);
    denom = 1 + (REV/Ki_bmal_rev)^hill_bmal_rev + (ROR/Ka_bmal_ror)^hill_bmal_ror;
    bmal_new = bmal + (-dm_bmal * bmal + numer/denom) * dt;
    
    % 6
    PER_new = PER + (-dp_per *(1+m_per_sirt*Act_SIRT + m_per_ampk*Act_AMPK)*PER + kp_per * per) * dt;

    % 7
    CRY_new = CRY + (-dp_cry*(1+m_cry_ampk*Act_AMPK)*CRY + kp_cry * cry ...
        - (kass_pc * CRY * PER - kdiss_pc * PC)) *dt;
    
    % 8
    REV_new = REV + (-dp_rev * REV + kp_rev * rev) * dt;

    % 9
    ROR_new = ROR + (-dp_ror * ROR + kp_ror * ror) * dt;

    % 10
    BMAL_new = BMAL + (-dp_bmal*BMAL + kp_bmal *bmal - (kass_cb * BMAL) - kdiss_cb * CB) * dt;

    % 11
    PC_new = PC + (kass_pc * CRY * PER - kdiss_pc * PC - d_pc * PC) * dt;

    % 12
    CB_new = CB + (kass_cb * BMAL - kdiss_cb * CB - d_cb * CB) * dt;

    % 13
    numer = Vmax_nampt * (1+fold_nampt * (CB/Ka_nampt_cb/(1+Act_SIRT))^hill_nampt_cb);
    denom = 1 + (CB/Ka_nampt_cb*(1+Act_SIRT))^hill_nampt_cb * (1 + (PC/Ki_nampt_pc)^hill_nampt_pc);
    nampt_new = nampt + (-dm_nampt * nampt + numer/denom) * dt;

    % 14
    NAMPT_new = NAMPT + (-(dp_nampt*NAMPT)/(1+m_nampt_ampk*Act_AMPK) + kp_nampt * nampt) * dt;

    % 15
    numer = Vnad* NAMPT * (NAD_tot-NAD);
    denom = Knam + NAD_tot - NAD;
    NAD_new = NAD + (-(d_nad*(NAD-NAD_basal)/(Knad+NAD-NAD_basal)) + numer/denom) * dt;

    % 16
    numer = Vmax_dbp * (1 + fold_dbp* (CB/Ka_dbp_cb/(1+Act_SIRT))^hill_dbp_cb);
    denom = 1 + (CB/Ka_dbp_cb/(1+Act_SIRT))^hill_dbp_cb * (1 + (PC/Ki_dbp_pc)^hill_dbp_pc);
    dbp_new = dbp + (numer/denom - dm_dbp*dbp)*dt;

    % 17
    Act_SIRT_new = Csirt*Vsirt*NAD/(Ksirt + NAD);
 
    % 18
    Act_AMPK_new = Campk*(amp1*ampk_vec(clock) + amp2 * ampk_vec(clock)) ...
        + (1-Campk) * offs;
    % 19
    numer = 1 * Vpg * Act_AMPK * Act_SIRT * PGC1a;
    denom = 1 + (Act_AMPK/Kpg1)*(1+ Act_SIRT/Kpg2);
    Act_PGC1a_new = numer / denom;

    % 20
    %PGC1a_new = amp3 * P(clock*dt, tc3, Td3, Tr3, 24);
    PGC1a_new = amp3 * ampk_vec(clock);
 
    % 21
    %agonist_rev_new = amp4 * P(clock*dt, tc4, Td4, Tr4, 24);
    agonist_rev_new = amp4 * ampk_vec(clock);

    % 22
    Ki_bmal_rev_new = (Ki_bmal_rev0) / (1 + agonist_rev);

    % 23
    Ki_cry_rev_new = (Ki_cry_rev0) / (1 + agonist_rev);
    
    % additionally
    
    if with_light == 1
        lighting_per_max_new = Vmax_per;% + Vmax_per * epsilon_vec(clock);
        lighting_cry_max_new = Vmax_cry;% + Vmax_cry * epsilon_vec(clock);
        lighting_ror_max_new = Vmax_ror + Vmax_ror * epsilon_vec(clock);
        lighting_rev_max_new = Vmax_rev + Vmax_rev * epsilon_vec(clock);
    elseif with_light == 0
        lighting_per_max_new = Vmax_per;% + Vmax_per * epsilon_vec(clock);
        lighting_cry_max_new = Vmax_cry;% + Vmax_cry * epsilon_vec(clock);
        lighting_ror_max_new = Vmax_ror;
        lighting_rev_max_new = Vmax_rev;
    end
    
        
    
    
    per= per_new; 
    cry= cry_new;        
    rev= rev_new;        
    ror= ror_new;         
    bmal= bmal_new;        
    PER= PER_new;
    CRY= CRY_new;
    REV= REV_new;
    ROR= ROR_new;
    BMAL= BMAL_new;
    PC= PC_new;
    CB= CB_new;
    nampt= nampt_new;
    NAMPT= NAMPT_new;
    dbp = dbp_new;
    NAD= NAD_new;
    Act_SIRT= Act_SIRT_new;
    Act_PGC1a=0;
    Act_AMPK = Act_AMPK_new;
    PGC1a= PGC1a_new;
    agonist_rev= agonist_rev_new;
    Ki_bmal_rev = Ki_bmal_rev_new;
    Ki_cry_rev = Ki_cry_rev_new;
    lighting_per_max = lighting_per_max_new;
    lighting_cry_max = lighting_per_max_new;
    lighting_rev_max = lighting_rev_max_new;
    lighting_ror_max = lighting_ror_max_new;
    
    
    per_vec(clock) = per_new;
    cry_vec(clock) = cry_new;
    rev_vec(clock) = rev_new;
    ror_vec(clock) = ror_new;
    bmal_vec(clock) = bmal_new;
    PER_vec(clock) = PER_new;
    CRY_vec(clock) = CRY_new;
    REV_vec(clock) = REV_new;
    ROR_vec(clock) = ROR_new;
    BMAL_vec(clock) = BMAL_new;
    PC_vec(clock) = PC_new;
    CB_vec(clock) = CB_new;
    nampt_vec(clock) = nampt_new;
    NAMPT_vec(clock) = NAMPT_new;
    dbp_vec(clock) = dbp_new;
    NAD_vec(clock) = NAD_new;
    Act_SIRT_vec(clock) = Act_SIRT_new;
    Act_PGC1a_vec(clock) = Act_PGC1a_new;
    Act_AMPK_vec(clock) = Act_AMPK_new;
    PGC1a_vec(clock) = PGC1a_new;
    agonist_rev_vec(clock) = agonist_rev_new;
    Ki_bmal_rev_vec(clock) = Ki_bmal_rev_new;
    Ki_cry_rev_vec(clock) = Ki_cry_rev_new;
    lighting_per_max_vec(clock) = lighting_per_max_new;
    lighting_cry_max_vec(clock) = lighting_cry_max_new;
    lighting_rev_max_vec(clock) = lighting_rev_max_new;
    lighting_ror_max_vec(clock) = lighting_ror_max_new;
    
end

%{
function p_value = P(t, t_c, T_d, T_w, T_c)
    % p oscillates between 0 and 1
    %p_value = S(t, t_c - T_d / 2, 0, T_w, T_c) - S(t, t_c - T_d /2, T_d, T_w, T_c) ...
    %    + S(t,t_c - T_d/2, T_c, T_w, T_c) - S(t, t_c - T_d/2, T_c + T_d, T_w, T_c);
    p_value = 1-meals(t);
    %p_value = ampk_vec(round(t/0.01));
    %p_value = 2; 
end


function s_value = S(t, t_f, t_s, T_w, T_c)
    s_value = (1 + tanh((tau(t-t_f, T_c)-t_s)/T_w))/2;
end

function tau_value = tau(t, T_c)
    tau_value = t - T_c * floor(t/T_c);
end
%}

end