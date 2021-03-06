addpath('./code/bottleneck')
addpath('./code/circuit_sim')

%% 1.
v=linspace(0.1,10,10);

%avg_currents=zeros(1,length(v));
%for k=1:length(v)
%    avg_currents(k)=coupled_sim("3B",v(k),0.2e-7);
%end

% Allow me to save you (and me) some time
avg_currents=[0.0330    0.3191    0.5672    0.7803    1.0331    1.2361    1.4916    1.7967    1.9779    2.0885];

figure();
plot(v,avg_currents);
title('Average Current vs. Voltage');
ylabel('Average Current (A)');
xlabel('Voltage (V)');   

disp("Done Q1. Press any key to continue...");
pause;
close all;

%% 2.
p=polyfit(v,avg_currents,1);
r3=1/p(1);

fprintf("R3: %e\n",r3);
disp("Done Q2. Press any key to continue...");
pause;

%% 3.
freq_domain_sim
disp("Done Q3. Press any key to continue...");
pause;

%% 4.
transient_sim
disp("Done Q4. Press any key to continue...");
pause;

%% 5.
transient_sim_II
disp("Done Q5. Press any key to continue...");
pause;

