clear;
close all;

[G,C]=circuitModel();
C1 = 0.25;

%% DC Case
V3 = [];
VO = [];
Vin=linspace(-10,10,20);
for v=Vin
    F = [v,0,0,0,0,0];
    V = G\F';
    V3(end+1) = V(3);
    VO(end+1) = V(5);
end

figure();
plot(Vin,V3,'-',Vin,VO,'-');
title('Voltage at Node vs. Input Voltage');
legend('V3','VO');
xlabel('Input Voltage (V)')
ylabel('Voltage at Node(V)')

%% AC Case
F = [1,0,0,0,0,0];
%afreq = linspace(0,10e9,1e5);
afreq = logspace(0,3,100);
VO = [];
gain = [];
for w = afreq
    V=(G+j*w*C)\F';
    VO(end+1) = V(5);
    gain(end+1) = 20*log(V(5));
end

figure();
subplot(2,1,1);
semilogx(afreq,VO);
title('Output Voltage vs. Frequency (rads/s)');
xlabel('Frequency (rads/s)')
ylabel('Output Voltage(V)')
subplot(2,1,2);
semilogx(afreq,gain);
title('Gain (dB) vs. Frequency (rads/s)');
xlabel('Frequency (rads/s)')
ylabel('Gain (dB)')

%% AC Case II
% Random changes to the capacitance
F = [1,0,0,0,0,0];
gain = zeros(1,500);
for i = 1:500
    C(2, 1) = normrnd(C1, 0.05); 
    C(2, 2) = normrnd(-C1, 0.05);    
    V=(G+j*pi*C)\F';    
    gain(i) = 20*log(V(5));
end

figure();
histogram(real(gain))
title('Gain (dB) for Normally Distributed Capacitance')
xlabel('Gain (dB)')
ylabel('# Occurrences')