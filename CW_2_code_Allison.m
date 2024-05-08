%% Question 2_1

%% Calculating the firing rate of coupled LIF Neurons with equal synaptic weights

clear all
close all

% First neuron specification

R_1 = 10 * 10^6;       % 10 Mega ohm 
C_1 = 0.1 * 10^-9;     % 0.1 Nano Farad
tau_1= R_1*C_1;
E_1 = -75 * 10^-3;     % resting potential 75 mili volt 
I_ext_1 = 2.501 * 10^-9;       % I is 2 nano amper
v_thr_1 = -50 * 10^-3;
v_spike_1 = 20 *10^-3;
w_1 = 2 * 10^-9;

% Second neuron specification

R_2 = 10 * 10^6;       % 10 Mega ohm 
C_2 = 0.1 * 10^-9;     % 0.1 Nano Farad
tau_2= R_2*C_2;
E_2 = -75 * 10^-3;     % resting potential 75 mili volt 
I_ext_2 = 2.501 * 10^-9;       % I is 2 nano amper
v_thr_2 = -50 * 10^-3;
v_spike_2 = 20 *10^-3;
w_2 = 2 * 10^-9;


% simulation time 
dt=0.1/1000;
t=[0:dt:100/1000] ; % t is mili second 

% response of presynaptic neuron 
if_spike_1 = zeros(size(t));

% response of presynaptic neuron 
if_spike_2 = zeros(size(t));

% numerical solution -> Euler Method 

v2(1)=E_2;
v1(1)=E_1;

 
for l=2:length(t)

% calculation of first neuron membrain potential     
I_syn_1 = w_1 * if_spike_2(l-1);
I_1 = I_ext_1 + I_syn_1 ;   

v1(l) = v1(l-1) + dt*((E_1-v1(l-1))/tau_1 + I_1/C_1);

if v1(l)>v_thr_1
    v1(l) = v_spike_1;
end 

if v1(l-1)==v_spike_1
    v1(l) = E_1;
    if_spike_1(l) = 1 ;
end



% calculation of second neuron membrain potential     
I_syn_2 = w_2 * if_spike_1(l-1);
I_2 = I_ext_2 + I_syn_2 ;   

v2(l) = v2(l-1) + dt*((E_2-v2(l-1))/tau_2 + I_2/C_2);

if v2(l)>v_thr_2
    v2(l) = v_spike_2;
end 

if v2(l-1)==v_spike_2
    v2(l) = E_2;
    if_spike_2(l) = 1; 
end


end


%Firing rates

out_spikes_1 = (v1 == v_spike_1);
neuron1_fr = sum(out_spikes_1)/t(end);

out_spikes_2 = (v2 == v_spike_2);
neuron2_fr = sum(out_spikes_2)/t(end);


%figure(1)
%plot(t,v1)
%hold on
%plot(t,v2)

%xlabel("Time (s)")
%ylabel("Membrane potential of N1/N2 (v)")
%title("LIF Synaptic input - Allison 33792133")



%% Calculating the firing rate of coupled LIF Neurons with different synaptic weights


% First neuron specification

R_3 = 10 * 10^6;       % 10 Mega ohm 
C_3 = 0.1 * 10^-9;     % 0.1 Nano Farad
tau_3= R_3*C_3;
E_3 = -75 * 10^-3;     % resting potential 75 mili volt 
I_ext_3 = 2.501 * 10^-9;       % I is 2 nano amper
v_thr_3 = -50 * 10^-3;
v_spike_3 = 20 *10^-3;
w_3 = 4 * 10^-9;

% Second neuron specification

R_4 = 10 * 10^6;       % 10 Mega ohm 
C_4 = 0.1 * 10^-9;     % 0.1 Nano Farad
tau_4= R_4*C_4;
E_4 = -75 * 10^-3;     % resting potential 75 mili volt 
I_ext_4 = 2.501 * 10^-9;       % I is 2 nano amper
v_thr_4 = -50 * 10^-3;
v_spike_4 = 20 *10^-3;
w_4 = 1 * 10^-9;


% simulation time 
dt=0.1/1000;
t=[0:dt:100/1000] ; % t is mili second 

% response of presynaptic neuron 
if_spike_3 = zeros(size(t));

% response of presynaptic neuron 
if_spike_4 = zeros(size(t));

% numerical solution -> Euler Method 

v4(1)=E_4;
v3(1)=E_3;

 
for l=2:length(t)

% calculation of first neuron membrain potential     
I_syn_3 = w_3 * if_spike_4(l-1);
I_3 = I_ext_3 + I_syn_3 ;   

v3(l) = v3(l-1) + dt*((E_3-v3(l-1))/tau_3 + I_3/C_3);

if v3(l)>v_thr_3
    v3(l) = v_spike_3;
end 

if v3(l-1)==v_spike_3
    v3(l) = E_3;
    if_spike_3(l) = 1 ;
end



% calculation of second neuron membrain potential     
I_syn_4 = w_4 * if_spike_3(l-1);
I_4 = I_ext_4 + I_syn_4 ;   

v4(l) = v4(l-1) + dt*((E_4-v4(l-1))/tau_4 + I_4/C_4);

if v4(l)>v_thr_4
    v4(l) = v_spike_4;
end 

if v4(l-1)==v_spike_4
    v4(l) = E_4;
    if_spike_4(l) = 1; 

end


end


%Firing rates

out_spikes_3 = (v3 == v_spike_3);
neuron3_fr = sum(out_spikes_3)/t(end);


out_spikes_4 = (v4 == v_spike_4);
neuron4_fr = sum(out_spikes_4)/t(end);



%figure(1)
%plot(t,v1)
%hold on
%plot(t,v2)

%xlabel("Time (s)")
%ylabel("Membrane potential of N1/N2 (v)")
%title("LIF Synaptic input - Allison 33792133")


% printing resulting firing rates, question 2_1 answer
fprintf('\n\nQuestion 2_1 Printed Output \n\n')

fprintf('Equal Synaptic Weights- \n')
fprintf('Neuron 1 (w_1 = 2): %i' , neuron1_fr)
fprintf(' Hz\nNeuron 2 (w_2 = 2): %i' , neuron2_fr)
fprintf(' Hz\n\n')

fprintf('Different Synaptic Weights- \n')
fprintf('Neuron 3 (w_1 = 4): %i' , neuron3_fr)
fprintf(' Hz\nNeuron 4 (w_2 = 1): %i' , neuron4_fr)
fprintf(' Hz\n')



%% Question 2_2


clear all
close all

%neuron 1

R_1 = 10 * 10^6;     
C_1 = 0.1 * 10^-9;    
E_1 = -75 * 10^-3;     
I_ext_1 = 2.501 * 10^-9;     
v_thr_1 = -50 * 10^-3;
v_spike_1 = 20 *10^-3;
tau_1= R_1*C_1;

%neuron 2

R_2 = 10 * 10^6;     
C_2 = 0.1 * 10^-9;    
E_2 = -75 * 10^-3;     
I_ext_2 = 2.501 * 10^-9;      
v_thr_2 = -50 * 10^-3;
v_spike_2 = 20 *10^-3;
tau_2= R_2*C_2;

%neuron 3

R_3 = 10 * 10^6;     
C_3 = 1 * 10^-9;    
E_3 = -75 * 10^-3;     
I_ext_3 = 2.501 * 10^-9;      
v_thr_3 = -50 * 10^-3;
v_spike_3 = 20 *10^-3;
tau_3= R_3*C_3;


%synaptic weights

w_1_2 = 0.1 * 10^-9;
w_2_1 = 0.1 * 10^-9;
w_1_3 = 50 * 10^-9;
w_2_3 = 0.2 * 10^-9;


% simulation time 
dt=0.1/1000;
t=[0:dt:100/1000] ; % t is mili second 

% response of presynaptic neuron 
if_spike_1 = zeros(size(t));

% response of presynaptic neuron 
if_spike_2 = zeros(size(t));

% response of presynaptic neuron 
if_spike_3 = zeros(size(t));


% numerical solution -> Euler Method 

v1(1)=E_1;
v2(1)=E_2;
v3(1)=E_3;


for l=2:length(t)

% calculation of first neuron membrain potential     
I_syn_1 = w_2_1 * if_spike_2(l-1);
I_1 = I_ext_1 + I_syn_1 ;   

v1(l) = v1(l-1) + dt*((E_1-v1(l-1))/tau_1 + I_1/C_1);

if v1(l)>v_thr_1
    v1(l) = v_spike_1;
end 

if v1(l-1)==v_spike_1
    v1(l) = E_1;
    if_spike_1(l) = 1 ;
end



% calculation of second neuron membrain potential     
I_syn_2 = w_1_2 * if_spike_1(l-1);
I_2 = I_ext_2 + I_syn_2 ;   

v2(l) = v2(l-1) + dt*((E_2-v2(l-1))/tau_2 + I_2/C_2);

if v2(l)>v_thr_2
    v2(l) = v_spike_2;
end 

if v2(l-1)==v_spike_2
    v2(l) = E_2;
    if_spike_2(l) = 1; 
end


% calculation of third neuron membrain potential     
I_syn_3 = w_1_3 * if_spike_1(l-1) + w_2_3 * if_spike_2(l-1);
I_3 = I_ext_3 + I_syn_3 ;   

v3(l) = v3(l-1) + dt*((E_3-v3(l-1))/tau_3 + I_3/C_3);

if v3(l)>v_thr_3
    v3(l) = v_spike_3;
end 

if v3(l-1)==v_spike_3
    v3(l) = E_3;
    if_spike_3(l) = 1; 
end


end


figure(1)
plot(t,v1)
hold on
plot(t,v2)
hold on 
plot(t,v3)

xlabel("Time (s)")
ylabel("Membrane potential of N1/N2/N3 (v)")
title("LIF Synaptic input three coupled neurons- Allison")


