% Define parameters
Gl = 30e-9;  
Gnamax = 12e-6;
Gkmax = 3.6e-6;
Ena = 45e-3;
Ek = -82e-3;
El = -60e-3;
Cm = 100e-12;
Iapp = 0;

% Define time parameters
t_start = 0;  % Start time
t_end = 0.35; % End time
dt = 0.00001;   % Time step size

% Define time vector
t = t_start:dt:t_end;

% Preallocate arrays to store results
V = zeros(size(t));
m = zeros(size(t));
h = zeros(size(t));
n = zeros(size(t));

% Initial conditions
V(1) = El;

% Main simulation loop
% Matlab starts at index 1
% End at length(t)-1 because you've got i+1 indexes
for i = 1:length(t)-1
    if V(i)==-0.045
        am=100;
    else
        am=(1e5*(-V(i)-0.045))/(exp(100*(-V(i)-0.045))-1);
    end
    bm=4e3*exp((-V(i)-0.070)/(0.018));
    ah=70*exp(50*(-V(i)-0.070));
    bh=1e3/(1+exp(100*(-V(i)-0.040)));
    if V(i)==-0.06
        an=100;
    else
        an=(1e4*(-V(i)-0.060))/(exp(100*(-V(i)-0.060))-1);
    end
    bn=125*exp((-V(i)-0.070)/(0.08));
    dV = (Gl*(El-V(i))+Gnamax*power(m(i),3)*h(i)*(Ena-V(i))+Gkmax*power(n(i),4)*(Ek-V(i))+Iapp)/Cm;
    dm = am*(1-m(i))- bm*m(i);
    m(i+1) = m(i) + dm * dt;
    dh = ah*(1-h(i))-bh*h(i);
    h(i+1) = h(i) + dh * dt;
    dn = an*(1-n(i))-bn*n(i);
    n(i+1) = n(i) + dn * dt;
    V(i+1) = V(i) + dt*dV;
end

% Plot results
figure(1)
plot(t, V);
xlabel('Simulation Time (s)');
ylabel('Membrane Potential (V)');
title('Membrane Potential vs. Time');
