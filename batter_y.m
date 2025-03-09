%% Renewable source integration to grid
% PV array one stage DC/AC VSC interface
% 
% Created: 15.11.2024
% Version: VSC AC Side connected to ideal grid
% Control structure: Inner Current Control Loop
%                    Phase-Locked Loop
%                    I/V abc to qd0 transformation
%                    id (reactive power reference)
%                    iq (active power reference)
%
% Mod: 27.11.2024
% Version: PV Array + DC/AC VSC connected to ideal grid
% Control structure: MPPT (open circuit method)
%                    Outer DC Voltage Loop
%                    Outer Reactive Power Loop
%                    Inner Current Control Loop
%                    Phase-Locked Loop
%                    I/V abc to qd0 transformation
%                    id (reactive power reference)
%                    iq (active power reference)
%
% Mod: 03.12.2024
% Included P-V curves of PV array
%
% Author: Ullón Huber, Mauricio









% Simulation data for the IEEE 13 Node Test Feeder model

% miles/km
mi2km = 1.60934;

% feet to km
ft2km = 0.0003048;

% microsiemens to Farads
ms2F = 1/2/pi/60*1e-6;


%% Configuration 601 - series reactance - ohm/mile

R_601 = [0.3465 0.1560 0.1580;0.1560 0.3375 0.1535;0.1580 0.1535 0.3414];
X_601 = [1.0179 0.5017 0.4236;0.5017 1.0478 0.3849;0.4236 0.3849 1.0348];

% charging susceptance - microsiemens/mile

B_601 = [6.2998 -1.9958 -1.2595;-1.9958 5.9597 -0.7417;-1.2595 -0.7417 5.6386];

% convert for SPS

R_601 = R_601/mi2km;

L_601 = X_601/mi2km/2/pi/60;

C_601 = B_601/mi2km*ms2F;


%% Configuration 602 - series reactance - ohm/mile

R_602 = [0.7526 0.1580 0.1560;0.1580 0.7475 0.1535;0.1560 0.1535 0.7436];
X_602 = [1.1814 0.4236 0.5017;0.4236 1.1983 0.3849;0.5017 0.3849 1.2112];

% charging susceptance - microsiemens/mile

B_602 = [5.6990 -1.0817 -1.6905;-1.0817 5.1795 -0.6588;-1.6905 -0.6588 5.4246];

% convert for SPS

R_602 = R_602/mi2km;

L_602 = X_602/mi2km/2/pi/60;

C_602 = B_602/mi2km*ms2F;


%% Configuration 603 - series reactance - ohm/mile

R_603 = [0 0 0;0 1.3294 0.2066;0 0.2066 1.3238];
X_603 = [0 0 0;0 1.3471 0.4591;0 0.4591 1.3569];

% charging susceptance - microsiemens/mile

B_603 = [0 0 0;0 4.7097 -0.8999;0 -0.8999 4.6658];

% convert for SPS

R_603 = R_603/mi2km;

L_603 = X_603/mi2km/2/pi/60;

C_603 = B_603/mi2km*ms2F;


%% Configuration 604 - series reactance - ohm/mile

R_604 = [1.3238 0 0.2066;0 0 0;0.2066 0 1.3294];
X_604 = [1.3569 0 0.4591;0 0 0;0.4591 0 1.3471];

% charging susceptance - microsiemens/mile

B_604 = [4.6658 0 -0.8999;0 0 0;-0.8999 0 4.7097];


% convert for SPS

R_604 = R_604/mi2km;

L_604 = X_604/mi2km/2/pi/60;

C_604 = B_604/mi2km*ms2F;


%% Configuration 605 - series reactance - ohm/mile

R_605 = [0 0 0;0 0 0;0 0 1.3292];
X_605 = [0 0 0;0 0 0;0 0 1.3475];

% charging susceptance - microsiemens/mile

B_605 = [0 0 0;0 0 0;0 0 4.5193];


% convert for SPS

R_605 = R_605/mi2km;

L_605 = X_605/mi2km/2/pi/60;

C_605 = B_605/mi2km*ms2F;

%% Configuration 606 - series reactance - ohm/mile

R_606 = [0.7982 0.3192 0.2849;0.3192 0.7891 0.3192;0.2849 0.3192 0.7982];
X_606 = [0.4463 0.0328 0.0143;0.0328 0.4041 0.0328;0.0143 0.0328 0.4463];
%X_606 = [0.4463 0.0328 -0.0143;0.0328 0.4041 0.0328;-0.0143 0.0328 0.4463];

% charging susceptance - microsiemens/mile

B_606 = [96.8897 -1e-6 -1e-6;-1e-6 96.8897 -1e-6;-1e-6 -1e-6 96.8897];

% convert for SPS

R_606 = R_606/mi2km;

L_606 = X_606/mi2km/2/pi/60;

C_606 = B_606/mi2km*ms2F;


%% Configuration 607 - series reactance - ohm/mile

R_607 = [1.3425 0 0;0 0 0;0 0 0];
X_607 = [0.5124 0 0;0 0 0;0 0 0];

% charging susceptance - microsiemens/mile

B_607 = [88.9912 0 0;0 0 0;0 0 0];

% convert for SPS

R_607 = R_607/mi2km;

L_607 = X_607/mi2km/2/pi/60;

C_607 = B_607/mi2km*ms2F;

%% PV module (Model: KC200GT Solar)
Imp_n=7.61;                 % Nominal max. power current
Vmp_n=26.3;                 % Nominal max. power voltage
Pmax_e=200.143;             % Nominal max. power
Isc_n=8.21;                 % Nominal short-circuit current
Voc_n=32.9;                 % Nominal open-circuit voltage
KV=-0.123;                  % Temp. coeff. of Voc
KVmp=Vmp_n/Voc_n;           % Ratio max. power to open-circuit (MPPT)
KI=0.0032;                  % Temp. coeff. of Isc
Ns_cell=54;                 % Number of cells connected in series
a=1.3;                      % Diode ideality constant
Rp=425.405;                 % Eq. parallel resistance
Rs=0.221;                   % Eq. series resistance

% Size of the PV Array
Ns=45;                      % Number of modules in series
Np=345;                     % Number of modules in parallel

% Standard test conditions (STC)
Gn=1000;                    % Irradiance
Tn=25;                      % Temperature

% Aditional values for calculations
k=1.3806503e-23;            % Bolzmann's constant
q=1.6021764e-19;            % Charge of electron

% Nominal Ipv considering Rp and Rs
Ipv_n=Isc_n*(Rp+Rs)/Rp;



%% Grid parameters
Ug=480;                     % Grid voltage (L2L-rms)
Vg=Ug/sqrt(3);              % L2L-rms to L2N-rms
fg=60;                      % Nominal frequency
Vp=Ug*sqrt(2)/sqrt(3);      % Peak voltage (L2G)

%% VSC parameters
Scn=3e6;                    % Nominal converter power
Ucn=480;                    % Nominal converter voltage (L2L-rms)
Vdc=1200;                   % Nominal DC voltage
Cdc=278e-3;                 % Capacitance of VSC DC bus
Rc=0.768e-3;                % Resistance of coupling filter
Lc=20.372e-6;               % Inductance of coupling filter
ma=1;                       % Modulation index

% Converter limits
Ic_max=Scn/(3*Vg);              % Current
Vc_max=(1/sqrt(2))*ma*(Vdc/2);  % Voltage (L2N-rms)
Uc_max=Vc_max*sqrt(3);          % Voltage (L2L-rms)

fprintf('VSC output limits\n');
fprintf('Current: Ic_max = %.2f A\n',Ic_max);
fprintf('Voltage: Uc_max = %.2f V\n',Uc_max);

%% Grid impedance
% XLratio=10;                 % Ratio between reactance and resistance X/R
% SCR=5;                      % Short-circuit ratio
% Sscg=SCR*Scn;               % Short-circuit power of the grid
% Zsc=Ug^2/Sscg;              % Short-circuit impedance
% Rg=Zsc/sqrt(1+XLratio^2);   % Resistance of the grid
% Xg=XLratio*Rg;              % Reactance of the grid
% Lg=Xg/(2*pi*fg);            % Inductance of the grid
              
%% VSC Control parameters

% PLL
xi_PLL=sqrt(2)/2;
omega_PLL=2*pi*fg;
Kp_PLL=xi_PLL*2*omega_PLL/Vp;
tau_PLL=2*xi_PLL/omega_PLL;
Ki_PLL=Kp_PLL/tau_PLL;

% Current control loop
tau_C=1e-3;
Kp_C=Lc/tau_C;
Ki_C=Rc/tau_C;

% Power control loop
Ki_pq=3*Vp/2;
tau_P=10*tau_C;
Kp_P=tau_C/(Ki_pq*tau_P);
Ki_P=1/(Ki_pq*tau_P);

% DC voltage Control
xi_DC=sqrt(2)/2;
omega_DC=4/(0.1*xi_DC);
Kp_DC=2*xi_DC*omega_DC*Cdc;
Ki_DC=omega_DC^2*Cdc;
Ts=5e-5;

%% Initialization Section
% System and Inverter Parameters
Pmod = 15;                     % Number of modules
Pnom = Pmod * 100e3;           % Nominal power (W)
Vnom_dc = 480;                 % Nominal DC link voltage (V)
Vnom_prim = 2400;              % Nominal inverter primary line-to-line voltage (Vrms)
Vnom_sec = 0.85 * Vnom_dc / 2 / sqrt(2) * sqrt(3); % Nominal inverter secondary line-to-line voltage (Vrms)
Fnom = 60;                     % Nominal frequency (Hz)
Bcap = 500;                    % Backup capacity (units)
TotalLeakage = 0.06;           % Transformer total leakage (pu)

% Transformer Parameters
Rxf = TotalLeakage / 25;       % Transformer resistance (pu)
Lxf = TotalLeakage;            % Transformer inductance (pu)

% Control Parameters
Ts_Control = 5e-5;             % Control sampling time (s)

% Inverter Choke RL Parameters
Rchoke = 0.15 / 100;           % Choke resistance (pu)
Lchoke = 0.15;                 % Choke inductance (pu)

% Base Impedance Calculation
Zbase = Vnom_sec^2 / Pnom;     % Base impedance (Ohm)

%% Filter Parameters
Qc = 0.1 * Pnom;               % Capacitive reactive power (VAR)
Pc = Qc / 50;                  % Active power (W)

% DC Link Capacitance
Clink = 2 * (3 / 4) * (Pnom / Fnom * 2 / Vnom_dc^2); % DC link energy for 3/4 cycle of Pnom

%% VDC Regulator Parameters
Kp_VDCreg = 2;                 % Proportional gain
Ki_VDCreg = 400;               % Integral gain

%% Feedforward Values
Rtot_pu = Rxf + Rchoke;        % Total resistance (pu)
Ltot_pu = Lxf + Lchoke;        % Total inductance (pu)

%% Current Regulator Parameters
Kp_Ireg = 0.3;                 % Proportional gain
Ki_Ireg = 20;                  % Integral gain


%% WIND TURBINE PARAMETERS


% Define constants and parameters
c1 = 0.005; c2 = 1.53; c3 = 0.5; c4 = 0.18; c5 = 121;
c6 = 27.9; c7 = 198; c8 = 2.36; c9 = 5.74; c10 = 11.35;
c11 = 16.1; c12 = 201;
rho = 1.225; % Air density (kg/m^3)
R = 63; % Rotor radius (m)
J = para_mdl.Jr + para_mdl.ngear^2 * para_mdl.Jg; % Combined inertia

% Wind turbine pitch angles (beta) and wind speeds (v)
beta = [0.0000, 0.0306, 0.0560, 0.0772, 0.0953, 0.1113, 0.1254, 0.1382, ...
        0.1499, 0.1607, 0.1783, 0.2134, 0.2421, 0.2669, 0.2890, 0.3091, ...
        0.3279, 0.3456, 0.3625, 0.3789]';
v = [11.2600, 11.5378, 11.8156, 12.0933, 12.3711, 12.6489, 12.9267, ...
     13.2044, 13.4822, 13.7600, 14.2600, 15.4533, 16.6467, 17.8400, ...
     19.0333, 20.2267, 21.4200, 22.6133, 23.8067, 25.0000]';
omega_r = 1.2671; % Rated rotor speed (rad/s)

% Define auxiliary functions
lambda_i = @(lambda, beta) 1 ./ (lambda + 0.08 * beta) - 0.035 ./ (c11 + c12 * beta.^3);
f1 = @(lambda) c4 ./ lambda;
f2 = @(lambda, beta) c5 * lambda_i(lambda, beta) - c6 * beta - c7 * beta.^c8 - c9;
f3 = @(lambda, beta) exp(-c10 * lambda_i(lambda, beta));

% Define derivatives of lambda_i
dlambda_i_dlambda = @(lambda, beta) -1 ./ (lambda + 0.08 * beta).^2;
dlambda_i_dbeta = @(lambda, beta) -0.08 ./ (lambda + 0.08 * beta).^2 + ...
                                   (3 * 0.035 * c12 * beta.^2) ./ (c11 + c12 * beta.^3).^2;

% Define analytical CQ function and its derivatives
cQ_tilde = @(lambda, beta) c1 * (1 + c2 * (beta + c3).^0.5) + ...
                          f1(lambda) .* f2(lambda, beta) .* f3(lambda, beta);
cQ = @(lambda, beta) cQ_tilde(lambda, beta) .* (1 + sign(cQ_tilde(lambda, beta))) / 2;

dcQ_dlambda = @(lambda, beta) -c4 ./ lambda.^2 .* f2(lambda, beta) .* f3(lambda, beta) + ...
                              f1(lambda) .* (c5 * dlambda_i_dlambda(lambda, beta)) .* f3(lambda, beta) + ...
                              f1(lambda) .* f2(lambda, beta) .* ...
                              (-c10 * exp(-c10 * lambda_i(lambda, beta)) .* dlambda_i_dlambda(lambda, beta));

dcQ_dbeta = @(lambda, beta) f1(lambda) .* (c5 * dlambda_i_dbeta(lambda, beta) - c6 - ...
                  c7 * c8 * beta.^(c8 - 1)) .* f3(lambda, beta) + ...
                  f1(lambda) .* f2(lambda, beta) .* ...
                  (-c10 * exp(-c10 * lambda_i(lambda, beta)) .* dlambda_i_dbeta(lambda, beta)) + ...
                  0.5 * c1 * c2 * (beta + c3).^(-0.5);

% Lambda as a function of wind speed and rotor speed
lambda = @(V, omega_r) omega_r * R ./ V;

% Initialize vectors for results
k_omega = zeros(length(v), 1);
k_beta = zeros(length(v), 1);
k_V = zeros(length(v), 1);

% Compute coefficients for each wind speed
for i = 1:length(v)
    k_omega(i) = 0.5 * rho * v(i)^2 * pi * R^3 * ...
                 (R / v(i) * dcQ_dlambda(lambda(v(i), omega_r), beta(i)));
    k_beta(i) = 0.5 * rho * v(i)^2 * pi * R^3 * ...
                dcQ_dbeta(lambda(v(i), omega_r), beta(i));
    k_V(i) = 0.5 * rho * pi * R^3 * ...
             (2 * v(i) * cQ(lambda(v(i), omega_r), beta(i)) - ...
              omega_r * R * dcQ_dlambda(lambda(v(i), omega_r), beta(i)));
end

% Compute control parameters
tau_ref = 4;
a = -k_omega / J;
b = k_beta / J;
k_P = 5./(b.*tau_ref);
k_I = a./(b.*tau_ref);

% Create a table with results and export to Excel
resultsTable = table(v, k_omega, k_beta, k_V, k_P, k_I, ...
    'VariableNames', {'WindSpeed', 'k_omega', 'k_beta', 'k_V', 'k_P', 'k_I'});
writetable(resultsTable, 'output1.xlsx');

