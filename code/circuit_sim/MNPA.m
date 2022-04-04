clear;
close all;

R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
Ro = 1000;
C1 = 0.25;
L = 0.2;
a = 100;

%% The G-matrix represents the voltage relationships bewteen the linear components of the circuit;
% It gets filled with the coeffcients of the set of equations 

% For each equation, let's write out the coeffcients in the following order: V1, V2, V3, V4, VO, IL
G=zeros(6);

% Node 1:
G(1,1)=1;

% Node 2:
G(2,1)=1/R1;
G(2,2)=-(1/R1+1/R2);
G(2,6)=-1;

% Node 3:
G(3,3)=1/R3;
G(3,6)=-1;

% Node 4:
G(4,3)=-a*1/R3;
G(4,4)=1;

% Node 5:
G(5,4)=1/R4;
G(5,5)=-(1/R4+1/Ro);

% The inductor:
G(6,2)=-1;
G(6,3)=1;

% The C-matrix represents the 1st order voltage relationships between the
% circuit components (I think)
C=zeros(6);

%From node 2:
C(2,1)=-C1;
C(2,2)=C1;

%From the inductor:
C(6,6)=L;

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

%% AC Case
F = [1,0,0,0,0,0];
afreq = linspace(0,100,50);
VO = [];
gain = [];
for w = afreq
    V=(G+j*w*C)\F';
    VO(end+1) = V(5);
    gain(end+1) = 20*log(V(5));
end

figure();
subplot(2,1,1);
plot(afreq,VO);
title('Output Voltage vs. Frequency (rads/s)');
subplot(2,1,2);
plot(afreq,gain);
title('Gain (dB) vs. Frequency (rads/s)');

%% AC Case II
% Random changes to the capacitance
F = [1,0,0,0,0,0];
gain = zeros(1,500);
for i = 1:500
    C(2, 1) = normrnd(-C1, 0.05); 
    C(2, 2) = normrnd(C1, 0.05);    
    V=(G+j*pi*C)\F';    
    gain(i) = 20*log(V(5));
end

figure();
histogram(real(gain))
title('Gain (dB) for Normally Distributed Capacitance')
