clear;
close all;

[G,C]=circuitModel_II(0.00001);

duration=1; %1 second simulation
num_steps=1000;
dt=duration/num_steps;

%% Gaussian pulse

t=(1:num_steps)*dt;
sd=0.03;
mu=5*sd; %3 sd to make the pulse start approximately at 0, 2 sd more to delay the pulse by 0.06 s
V_in=exp(-0.5*((t-mu)/sd).^2);

sim(num_steps,dt,G,C,V_in);


%% Again with 3 different capacitances

% smaller
[G,C]=circuitModel_II(0.00000001);
sim(num_steps,dt,G,C,V_in);

%bigger
[G,C]=circuitModel_II(0.001);
sim(num_steps,dt,G,C,V_in);

%more bigger
[G,C]=circuitModel_II(0.01);
sim(num_steps,dt,G,C,V_in);

%% Again with 2 different Time steps

[G,C]=circuitModel_II(0.00001);

%bigger
num_steps=100;
dt=duration/num_steps;
t=(1:num_steps)*dt;
V_in=exp(-0.5*((t-mu)/sd).^2);

sim(num_steps,dt,G,C,V_in);

%smaller
num_steps=10000;
dt=duration/num_steps;
t=(1:num_steps)*dt;
V_in=exp(-0.5*((t-mu)/sd).^2);

sim(num_steps,dt,G,C,V_in);

function [] = sim(num_steps,dt,G,C,V_in) 
    V_out=zeros(1,num_steps);
    V_prev=[0; 0; 0; 0; 0; 0];
    a=100;
    for k=1:num_steps
        % I excluded the noise current in the matrix formulation, so
        % they'll go in the forcing vector
        I_n=randn()*0.001;
        F=[V_in(k); 0; I_n; a*I_n; 0; 0];
        
        V=(G+C/dt)\(C*V_prev/dt+F);
        V_out(k)=V(5);
        V_prev=V;
    end

    figure();
    subplot(1,2,1);
    plot((1:num_steps)*dt,V_in,'-',(1:num_steps)*dt,V_out,'-');       
    title('Voltage vs. Time');
    legend('Vin','Vo');
    ylabel('Voltage (V)')
    xlabel('Time (s)')
    
    subplot(1,2,2);     
    omega=linspace(0.5*-1/dt,0.5*1/dt,num_steps);
    hold on;
    plot(omega,20*log10(abs(fftshift(fft(V_in)))));
    plot(omega,20*log10(abs(fftshift(fft(V_out)))));   
    title('Gain vs. Frequency');
    legend('Vin','Vo');
    xlabel('freq (Hz)')
    ylabel('Gain (dB)')
end

