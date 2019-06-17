%{
function result = delta(t)
n=10000;
sigma_n = 1/n;
result = 1/sqrt(2*pi*sigma_n) * exp(-t^2 / 2 / sigma_n);
end
%}

%{
function sti = stimulus2(t, meal_time, meal_amplitude)

sti=0;

for k = 1 : length(meal_time)
    sti = sti + 4.1* meal_amplitude(k) * delta(t-meal_time(k));

end
end
%}

%{

%}

