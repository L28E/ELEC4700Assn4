clear;
close all;

[G,C]=circuitModel();

duration=1; %1 second simulation
num_steps=1000;
dt=duration/num_steps;

%% A. A step that transistions from 0 to 1 at 0.03s.

V_in=zeros(1,num_steps);
V_in(0.3/dt:end)=1;

sim(num_steps,dt,G,C,V_in);

%% B. A sin(2πf t) function with f = 1/(0.03) 1/s. Try a few other frequencies. Comment.
% As expected of a lowpass filter, amplitude of the output decreases when
% we increase the frequency.

t=(1:num_steps)*dt;
f=1/0.03;
V_in=sin(2*pi*f*t);

sim(num_steps,dt,G,C,V_in);

%% C. A guassian pulse with a magnitude of 1, std dev. of 0.03s and a delay of 0.06s.

t=(1:num_steps)*dt;
sd=0.03;
mu=5*sd; %3 sd to make the pulse start approximately at 0, 2 sd more to delay the pulse by 0.06 s
V_in=exp(-0.5*((t-mu)/sd).^2);

sim(num_steps,dt,G,C,V_in);

function [] = sim(num_steps,dt,G,C,V_in) 
    V_out=zeros(1,num_steps);
    V_prev=[0; 0; 0; 0; 0; 0];
    for k=1:num_steps
        F=[V_in(k); 0; 0; 0; 0; 0];
        V=(G+C/dt)\(C*V_prev/dt+F);
        V_out(k)=V(5);
        V_prev=V;
    end

    figure();
    subplot(1,2,1);
    plot((1:num_steps)*dt,V_in,'-',(1:num_steps)*dt,V_out,'-');       
    title('Voltage vs. Time');
    legend('Vin','Vo');
    xlabel('Voltage (V)')
    ylabel('Time (s)')
    subplot(1,2,2); 
end

