% m = 360; %kg
% b = 2200; %N*sec/m

r_joint = 0.025;
IM = 360;
IL = 22.4;% I_arm = 0.014 kg-m^2, r_pivot = 0.025 m, IL = I_arm/r_pivot 
bM = 2200;% b_arm = 0.1 Nm-s/rad
bL = 160;
k = 350000; %N/m

Np = 4.0;
lead = 0.003;
N = 2*3.14159*Np/lead;
i2T = 0.0276; %nm/a
eff = 0.9;
beta1 = N * i2T * eff;
ktau = i2T/r_joint;

% % linear parameters
% IM = gearRatio^2 * 1.06 * 10^(-5)/r_knee^2;%360;
% IL = 0.02914/r_knee^2;%200;%0.16kgm^2;
% bM = 2/r_knee^2;%2000;
% bL = 1/r_knee^2;%1250;
% k = 331/r_knee^2;%350000; %rad/nm, or N/m
% beta1 = 3.21/r_knee;

% natural_freq = sqrt(k/m);
% damping_ratio = b/(2*sqrt(m*k));
% damped_freq = natural_freq*sqrt(1-damping_ratio^2);
% damped_freq_hz = damped_freq/2/pi;

%% plant model
%P = tf(k*N*i2T*eff,[m b k]); %tfr function from harmonic drive torque to spring torque (high z)
