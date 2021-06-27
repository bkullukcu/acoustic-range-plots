clear all
clear
clc
%% Coefficients
D = 4.6906779762565158E-8; %Flexural Rigidity
E1 = 2.500 * 10^9 ; %[Pa] Polyimide Young's Modulus
E2 = 27.8 * 10^9; %[Pa] Silver Young's modulus
E3 = 2.450 * 10^9 ; %[Pa] Piezoelectric Young's Modulus 
c = 343 ; % [m/s] Speed of sound in air medium
c0 = 343 ; % [m/s] Speed of sound in air medium
hh3 = 2.20*10^-6 ; % [m] Top of Piezo layer
hh1 = 1.50*10^-6 ; % [m] Top of Structural Layer
hh2 = 1.70*10^-6 ; % [m] Top of Bottom electrode
t3 = 0.50*10^-6 ; % [m] Piezoelectric Layer Thickness
t2 = 1.50 * 10 ^-5 ;% [m] Structural Layer Thickness
t1 = 0.20*10^-6 ;%[m] Bottom electrode
z3 = t1+t2+(t3/2) ;%1.95 * 10^-6 ; %[m] Piezo layer central axis
z1 = t1/2 ;% 0.75 * 10^-6 ; %[m] Structural layer central axis
z2 = t1 + t2/2 ;%0.10 * 10^-6 ; %[m] Bottom electrode central axis
zp = (((t1*z1*E1)+(t2*z2*E2)+(t3*z3*E3))/((t1*E1)+(t2*E2)+(t3*E3))); % Location of the neutral plane
dist1 = hh1 - zp ; % PI Distance to Neutral Axis
dist2 = hh2 - zp ; % Bottom Electrode Distance to Neutral Axis
dist3 = hh3 - zp ; % PVDF Distance to Neutral Axis
hr = linspace_round(0,100,10); % Relative Humidity
lambda00 = sqrt(10.22); % Eigenvalue of fundamental 01 vibration mode
a = 130*10^-6 ; % Plate radius [m]
ro3= 1780 ; % [kg/m^3] PVDF Density
ro1 = 1420 ; % [kg/m^3] PI Density
ro2 = 2700 ; % [kg/m^3] Aluminum Density
mueff = (ro1*t1)+(ro2*t2)+(ro3*t3); % Effective Mass Per Unit Area
wn = ((lambda00/a)^2)*sqrt(D/mueff) ; % Natural Freq.
lambda = c/wn; % Wavelength
v1 = 0.340 ;% Structural Layer Poisson's Ratio
v2 = 0.290;% Bottom electrode Poisson's Ratio
v3 = 0.340 ;% Piezo Layer Poisson's Ratio
plate1 = E1/(1-v1^2); % Plate Modulus for PI
plate2 = E2/(1-v2^2); % Plate Modulus for Bottom Electrode
plate3 = E3/(1-v3^2); % Plate Modulus for PVDF
r = 1*10^-2 ; % [m] Distance from the baffled piston
zdist = z3 - zp; % Piezo Mid-Layer Distance to Neutral Axis
F = linspace_round(0,1000000); %Frequency Range
Z0 = 419; % [Rayl] Characteristic Impedance of Air
T = 293.15; %Atmospheric Temperature in Kelvin
T0 = 293.15; % Reference Atmospheric Temperature in Kelvin
T01 = 273.16; % Ambient Temperature in Kelvin
Psat = 10^(-6.8346*((T01/T)^1.261)+4.6151); % Saturation Temperature
Ps = 1;
Ps0 = 1;
h=hr*((Psat/Ps0)/(Ps/Ps0)); % Humidity
%% Plots
hold on
xlabel("Absorbtion Coefficient a, [dB/100m*atm]");
ylabel("Frequency/Pressure, [Hz/atm]");
title("Attentuation for Different Freq. and Humid.");
%% 0 Humidity

Fr00 = (1/Ps0)*(24+4.04*(10^4)*h(1)*((0.02+h(1))/(0.391+h(1))));
FrN0 = (1/Ps0)*((T0/T)^(1/2))*(9+280*h(1)*exp(-4.17*(((T0/T)^(1/3)-1)))); 
alpha0 = @(F)((Ps/Ps0)*F^2)*(((1.84*(10^-11)*(sqrt(T/T0))))+((T/T0)^-(5/2))*(0.01278*((exp(-2239.1/T))/(Fr00+((F^2)/Fr00))))+(0.1068*((exp(-3352/T))/(FrN0+((F^2)/FrN0))))) ;% 
B0 = arrayfun(alpha0,F);
plot(B0)
%% 10 Humidity
Fr010 = (1/Ps0)*(24+4.04*(10^4)*h(2)*((0.02+h(2))/(0.391+h(2))));
FrN10 = (1/Ps0)*((T0/T)^(1/2))*(9+280*h(2)*exp(-4.17*(((T0/T)^(1/3)-1)))); 
alpha1 = @(F)((Ps/Ps0)*F^2)*(((1.84*(10^-11)*(sqrt(T/T0))))+((T/T0)^-(5/2))*(0.01278*((exp(-2239.1/T))/(Fr010+((F^2)/Fr010))))+(0.1068*((exp(-3352/T))/(FrN10+((F^2)/FrN10))))); % 
B1 = arrayfun(alpha1,F);
plot(B1)

%% 20 Humidity
Fr020 = (1/Ps0)*(24+4.04*(10^4)*h(3)*((0.02+h(3))/(0.391+h(3))));
FrN20 = (1/Ps0)*((T0/T)^(1/2))*(9+280*h(3)*exp(-4.17*(((T0/T)^(1/3)-1)))); 
alpha2 = @(F)((Ps/Ps0)*F^2)*(((1.84*(10^-11)*(sqrt(T/T0))))+((T/T0)^-(5/2))*(0.01278*((exp(-2239.1/T))/(Fr020+((F^2)/Fr020))))+(0.1068*((exp(-3352/T))/(FrN20+((F^2)/FrN20))))); % 
B2 = arrayfun(alpha2,F);
plot(B2)

%% 30 Humidity
Fr030 = (1/Ps0)*(24+4.04*(10^4)*h(4)*((0.02+h(4))/(0.391+h(4))));
FrN30 = (1/Ps0)*((T0/T)^(1/2))*(9+280*h(4)*exp(-4.17*(((T0/T)^(1/3)-1)))); 
alpha3 = @(F)((Ps/Ps0)*F^2)*(((1.84*(10^-11)*(sqrt(T/T0))))+((T/T0)^-(5/2))*(0.01278*((exp(-2239.1/T))/(Fr030+((F^2)/Fr030))))+(0.1068*((exp(-3352/T))/(FrN30+((F^2)/FrN30))))); % 
B3 = arrayfun(alpha3,F);
plot(B3)

%% 40 Humidity
Fr040 = (1/Ps0)*(24+4.04*(10^4)*h(5)*((0.02+h(5))/(0.391+h(5))));
FrN40 = (1/Ps0)*((T0/T)^(1/2))*(9+280*h(5)*exp(-4.17*(((T0/T)^(1/3)-1)))); 
alpha4 = @(F)((Ps/Ps0)*F^2)*(((1.84*(10^-11)*(sqrt(T/T0))))+((T/T0)^-(5/2))*(0.01278*((exp(-2239.1/T))/(Fr040+((F^2)/Fr040))))+(0.1068*((exp(-3352/T))/(FrN40+((F^2)/FrN40))))); % 
B4 = arrayfun(alpha4,F);
plot(B4)

%% 50 Humidity
Fr050 = (1/Ps0)*(24+4.04*(10^4)*h(6)*((0.02+h(6))/(0.391+h(6))));
FrN50 = (1/Ps0)*((T0/T)^(1/2))*(9+280*h(6)*exp(-4.17*(((T0/T)^(1/3)-1)))); 
alpha5 = @(F)((Ps/Ps0)*F^2)*(((1.84*(10^-11)*(sqrt(T/T0))))+((T/T0)^-(5/2))*(0.01278*((exp(-2239.1/T))/(Fr050+((F^2)/Fr050))))+(0.1068*((exp(-3352/T))/(FrN50+((F^2)/FrN50))))); % 
B5 = arrayfun(alpha5,F);
plot(B5)

%% 60 Humidity
Fr060 = (1/Ps0)*(24+4.04*(10^4)*h(7)*((0.02+h(7))/(0.391+h(7))));
FrN60 = (1/Ps0)*((T0/T)^(1/2))*(9+280*h(7)*exp(-4.17*(((T0/T)^(1/3)-1)))); 
alpha6 = @(F)((Ps/Ps0)*F^2)*(((1.84*(10^-11)*(sqrt(T/T0))))+((T/T0)^-(5/2))*(0.01278*((exp(-2239.1/T))/(Fr060+((F^2)/Fr060))))+(0.1068*((exp(-3352/T))/(FrN60+((F^2)/FrN60))))); % 
B6 = arrayfun(alpha6,F);
plot(B6)

%% 70 Humidity
Fr070 = (1/Ps0)*(24+4.04*(10^4)*h(8)*((0.02+h(8))/(0.391+h(8))));
FrN70 = (1/Ps0)*((T0/T)^(1/2))*(9+280*h(8)*exp(-4.17*(((T0/T)^(1/3)-1)))); 
alpha7 = @(F)((Ps/Ps0)*F^2)*(((1.84*(10^-11)*(sqrt(T/T0))))+((T/T0)^-(5/2))*(0.01278*((exp(-2239.1/T))/(Fr070+((F^2)/Fr070))))+(0.1068*((exp(-3352/T))/(FrN70+((F^2)/FrN70))))); %
B7 = arrayfun(alpha7,F);
plot(B7)

%% 80 Humidity
Fr080 = (1/Ps0)*(24+4.04*(10^4)*h(9)*((0.02+h(9))/(0.391+h(9))));
FrN80 = (1/Ps0)*((T0/T)^(1/2))*(9+280*h(9)*exp(-4.17*(((T0/T)^(1/3)-1)))); 
alpha8 = @(F)((Ps/Ps0)*F^2)*(((1.84*(10^-11)*(sqrt(T/T0))))+((T/T0)^-(5/2))*(0.01278*((exp(-2239.1/T))/(Fr080+((F^2)/Fr080))))+(0.1068*((exp(-3352/T))/(FrN80+((F^2)/FrN80))))); % 
B8 = arrayfun(alpha8,F);
plot(B8)

%% 90 Humidity
Fr090 = (1/Ps0)*(24+4.04*(10^4)*h(10)*((0.02+h(10))/(0.391+h(10))));
FrN90 = (1/Ps0)*((T0/T)^(1/2))*(9+280*h(10)*exp(-4.17*(((T0/T)^(1/3)-1)))); 
alpha9 = @(F)((Ps/Ps0)*F^2)*(((1.84*(10^-11)*(sqrt(T/T0))))+((T/T0)^-(5/2))*(0.01278*((exp(-2239.1/T))/(Fr090+((F^2)/Fr090))))+(0.1068*((exp(-3352/T))/(FrN90+((F^2)/FrN90))))); % 
B9 = arrayfun(alpha9,F);
plot(B9)
%% 100 Humidity
Fr0100 = (1/Ps0)*(24+4.04*(10^4)*h(11)*((0.02+h(11))/(0.391+h(11))));
FrN100 = (1/Ps0)*((T0/T)^(1/2))*(9+280*h(11)*exp(-4.17*(((T0/T)^(1/3)-1)))); 
alpha10 = @(F)((Ps/Ps0)*F^2)*(((1.84*(10^-11)*(sqrt(T/T0))))+((T/T0)^-(5/2))*(0.01278*((exp(-2239.1/T))/(Fr0100+((F^2)/Fr0100))))+(0.1068*((exp(-3352/T))/(FrN100+((F^2)/FrN100))))); %
B10 = arrayfun(alpha10,F);
plot(B10)

legend("% 0 Humidity","% 10 Humidity","% 20 Humidity","% 30 Humidity","% 40 Humidity","% 50 Humidity","% 60 Humidity","% 70 Humidity","% 80 Humidity","% 90 Humidity","% 100 Humidity");
