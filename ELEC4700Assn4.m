addpath('./code/bottleneck')

%% 1.
v=linspace(0.1,10,10);

%avg_currents=zeros(1,length(v));
%for k=1:length(v)
%    avg_currents(k)=coupled_sim("3B",v(k),0.2e-7);
%end

% Allow me to save you (and also me) some time
avg_currents=[6.78587371e-09 6.25654368e-08 1.08680734e-07 1.58160602e-07 2.02856379e-07 2.56760490e-07 3.08744450e-07 3.63797926e-07 3.78306750e-07 4.18630588e-07];

figure();
plot(v,avg_currents);
title('Average Current vs. Voltage');
ylabel('Average Current (A)');
xlabel('Voltage (V)');   

disp("Done. Press any key to continue...");
pause;
%close all;

%% 2.
p=polyfit(v,avg_currents,1);
r3=1/p(1);

disp(fprintf("R3: %f",r3));
disp("Done. Press any key to continue...");
pause;
