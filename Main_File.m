clear all; close all; clc;
tic
global Vref ws PC Ki Yabs Yang m n PL PL0 QL a0 x0 V0 alpha beta H Xd Xdp Xdpp Xq Xqp Xqpp Td0p ...
       Td0pp Tq0p Tq0pp Rs Xls Dm KA TA KE TE KF TF Ax Bx TCH TSV RD MH bb1 bb2 aa1 aa2

m = 3; n = 9;
baseMVA = 100; ws = 2*pi*60; ws_vector = ws * ones(1, m);  w0 = ws_vector;

%==========Creating Ybus and solving load flow using MATPOWER==============
mpopt = mpoption('PF_ALG', 1, 'ENFORCE_Q_LIMS', 0, 'VERBOSE', 0, 'OUT_ALL', 0);
[RESULTS1, SUCCESS] = runpf('IEEE9Bus', mpopt);
mpc = IEEE9Bus;
bus1 = RESULTS1.bus; branch1 = RESULTS1.branch; gen1 = RESULTS1.gen;
Ybus1 = makeYbus(baseMVA, bus1, branch1); Ybus1 = full(Ybus1);
%============================End of Ybus===================================
Ybus = Ybus1; Y = Ybus; bus = bus1; branch = branch1; gen = gen1;
Yabs = abs(Y); Yang = angle(Y);

IC1 = bus(:, 8);
IC2 = bus(:, 9);
gen0 = zeros(n, 1); gen0(1:m) = gen(1:m, 2); genP = gen0; IC3 = genP;
IC3 = IC3 / baseMVA;
gen0 = zeros(n, 1); gen0(1:m) = gen(1:m, 3); genQ = gen0; genQ = genQ + bus(:, 6); IC4 = genQ;
IC4 = IC4 / baseMVA;
IC5 = bus(:, 3); IC5 = IC5 / baseMVA;
IC6 = bus(:, 4); IC6 = IC6 / baseMVA;
IC = [IC1 IC2 IC3 IC4 IC5 IC6];
PL = IC(:, 5); QL = IC(:, 6);
PL0 = PL; % Store original PL
PG = IC(:, 3); QG = IC(:, 4);
TH0 = IC(:, 2) * pi / 180; TH0 = TH0'; V0 = IC(:, 1); V0 = V0';
VG0 = V0(1:m); THG0 = TH0(1:m);

%==========================Dynamic Parameters=============================
MD = [ % machine data
    23.640                 6.4000                3.0100;            % 1 -  H
    0.1460                 0.8958                1.3125;            % 2 -  Xd
    0.0608                 0.1198                0.1813;            % 3 -  Xdp
    0.0489                 0.0881                0.1133;            % 4 -  Xdpp
    0.0969                 0.8645                1.2578;            % 5 -  Xq
    0.0969                 0.1969                0.2500;            % 6 -  Xqp
    0.0396                 0.0887                0.0833;            % 7 -  Xqpp
    8.9600                 6.0000                5.8900;            % 8 -  Tdop
    0.1150                 0.0337                0.0420;            % 9 -  Td0pp
    0.3100                 0.5350                0.6000;            % 10-  Tqop
    0.0330                 0.0780                0.1875;            % 11-  Tq0pp
    0.0041                 0.0026                0.0035;            % 12-  RS
    0.1200                 0.1020                0.0750;            % 13-  Xls
    0.1*(2*23.64)/ws       0.2*(2*6.4)/ws      0.3*(2*3.01)/ws      % 14-  Dm
    ];
ED = [ % excitation data
    20.000 * ones(1, m);       % 1- KA
    0.2000 * ones(1, m);       % 2- TA
    1.0000 * ones(1, m);       % 3- KE
    0.3140 * ones(1, m);       % 4- TE
    0.0630 * ones(1, m);       % 5- KF
    0.3500 * ones(1, m);       % 6- TF
    0.0039 * ones(1, m);       % 7- Ax
    1.5550 * ones(1, m);       % 8- Bx
    ];
TD = [ % turbine data
    0.10 * ones(1, m);         % 1- TCH
    0.05 * ones(1, m);         % 2- TSV
    0.05 * ones(1, m);         % 3- RD
    ];
H     = MD(1,:);
Xd    = MD(2,:);
Xdp   = MD(3,:);
Xdpp  = MD(4,:);
Xq    = MD(5,:);
Xqp   = MD(6,:);
Xqpp  = MD(7,:);
Td0p  = MD(8,:);
Td0pp = MD(9,:);
Tq0p  = MD(10,:);
Tq0pp = MD(11,:);
Rs    = MD(12,:);
Xls   = MD(13,:);
Dm    = MD(14,:); % Keep original damping coefficients
KA    = ED(1,:);
TA    = ED(2,:);
KE    = ED(3,:);
TE    = ED(4,:);
KF    = ED(5,:);
TF    = ED(6,:);
Ax    = ED(7,:);
Bx    = ED(8,:);
TCH   = TD(1,:);
TSV   = TD(2,:);
RD    = TD(3,:); % Governor droop
MH = 2*H./ws;

% Integral gain for AGC
Ki = [3, 3, 3]; % Adjust these values as needed

%=================Initial Values Calculation (State Variables)=============
Vphasor  =  VG0.*exp(1i*(THG0));
Iphasor  =  conj((PG(1:m)' + 1i*QG(1:m)')./Vphasor);
E0       =  Vphasor + (Rs + 1i*Xq).*Iphasor;
Em       =  abs(E0);
D0       =  angle(E0);
Id0      =  real(Iphasor.*exp(-1i*(D0 - pi/2)));
Iq0      =  imag(Iphasor.*exp(-1i*(D0 - pi/2)));
Edp0     = (Xq - Xqp).*Iq0;
Si2q0    = (Xls - Xq).*Iq0;
Eqp0     =  Rs.*Iq0 + Xdp.*Id0 + V0(1:m).*cos(D0 - TH0(1:m));
Si1d0    =  Eqp0 - (Xdp - Xls).*Id0;
Efd0     =  Eqp0 + (Xd - Xdp).*Id0;
TM0      = ((Xdpp - Xls)./(Xdp - Xls)).*Eqp0.*Iq0 + ((Xdp - Xdpp)./(Xdp - Xls)).*Si1d0.*Iq0 + ...
           ((Xqpp - Xls)./(Xqp - Xls)).*Edp0.*Id0 - ((Xqp - Xqpp)./(Xqp - Xls)).*Si2q0.*Id0 + ...
           (Xqpp - Xdpp).*Id0.*Iq0;
VR0      = (KE + Ax.*exp(Bx.*Efd0)).*Efd0;
RF0      =  (KF./TF).*Efd0;
Vref     =  V0(1:m) + VR0./KA;
PSV0     =  TM0;
PC       =  PSV0; % Initial power command

% Initial value for xi (integral of frequency deviation)
xi0 = zeros(1, m);

%------------------------------
x0(1:m)          = Eqp0;
x0(m+1:2*m)      = Si1d0;
x0(2*m+1:3*m)    = Edp0;
x0(3*m+1:4*m)    = Si2q0;
x0(4*m+1:5*m)    = D0;
x0(5*m+1:6*m)    = ws_vector;
x0(6*m+1:7*m)    = Efd0;
x0(7*m+1:8*m)    = RF0;
x0(8*m+1:9*m)    = VR0;
x0(9*m+1:10*m)   = TM0;
x0(10*m+1:11*m)  = PSV0;
x0(11*m+1:12*m)  = xi0; % AGC integrator initial values
x01 = x0;                      % Initial values of state variables
a0 = [Id0 Iq0 V0 TH0];         % Initial values of algebraic variables
x0 = [x01 a0];                 % Final initial values of DAEs

%===========================Load parameters================================
alpha = 2; beta = alpha;

%============================Solving Es====================================
[bb1, aa1] = ndgrid(1:m, 1:n); % Generator buses
[bb2, aa2] = ndgrid(m+1:n, 1:n); % Load buses
M1 = ones(1, 12*m); M2 = zeros(1, 2*(m + n)); MT = [M1 M2]; % Mass Matrix updated for xi
M0 = diag(MT);
st = 1e-1;    tsim = 60; % Simulation time
tspan = [0, tsim];
opt = odeset('Mass', M0, 'MaxStep', 0.05, 'RelTol', 1e-6, 'AbsTol', 1e-8);

% Run the simulation
[t, x] = ode15s(@DAEs, tspan, x0, opt);

toc

%==============================Post-Processing=============================
% Update state variable extraction indices due to new state variable xi
Eqp   = x(:, 1:m);
Si1d  = x(:, m+1:2*m);
Edp   = x(:, 2*m+1:3*m);
Si2q  = x(:, 3*m+1:4*m);
Delta = x(:, 4*m+1:5*m);
w     = x(:, 5*m+1:6*m);
Efd   = x(:, 6*m+1:7*m);
RF    = x(:, 7*m+1:8*m);
VR    = x(:, 8*m+1:9*m);
TM    = x(:, 9*m+1:10*m);
PSV   = x(:, 10*m+1:11*m);
xi    = x(:, 11*m+1:12*m); % AGC integrator
Id    = x(:, 12*m+1:13*m);
Iq    = x(:, 13*m+1:14*m);
V     = x(:, 14*m+1:14*m+n);
TH    = x(:, 14*m+n+1:14*m+2*n);

%========================Compute Generator Outputs=========================
Pgen = zeros(length(t), m);
Qgen = zeros(length(t), m);
for k = 1:length(t)
    Vphasor = V(k, 1:m) .* exp(1i * TH(k, 1:m));
    Delta_k = Delta(k, :);
    Id_k = Id(k, :);
    Iq_k = Iq(k, :);
    Iphasor = (Id_k + 1i * Iq_k) .* exp(1i * (Delta_k - pi/2));
    Sgen_k = Vphasor .* conj(Iphasor);
    Pgen(k, :) = real(Sgen_k);
    Qgen(k, :) = imag(Sgen_k);
end

%========================Compute Load Over Time============================
PLt_over_time = zeros(length(t), n);
for k = 1:length(t)
    if t(k) >= 5 && t(k) <= 5.5
        PLt_over_time(k, :) = PL0' + (t(k) - 5)*(1.05*PL0' - PL0')/0.5;
    elseif t(k) > 5.5
        PLt_over_time(k, :) = 1.05 * PL0';
    else
        PLt_over_time(k, :) = PL0';
    end
end

%========================Compute Total Losses Over Time====================
Pload_over_time = sum(PLt_over_time, 2);
TotalLosses = sum(Pgen, 2) - Pload_over_time;

%========================Compute Line Flows and Losses=====================
branch = mpc.branch;
from_bus_indices = branch(:, 1);
to_bus_indices   = branch(:, 2);
r_branch = branch(:, 3);
x_branch = branch(:, 4);
Ybranch = 1 ./ (r_branch + 1i * x_branch);

num_branches = size(branch, 1);
LineFlows = zeros(length(t), num_branches);
BranchLosses = zeros(length(t), num_branches);

for k = 1:length(t)
    V_k = V(k, :) .* exp(1i * TH(k, :));
    V_from_k = V_k(from_bus_indices)';
    V_to_k   = V_k(to_bus_indices)';
    Ibranch_k = Ybranch .* (V_from_k - V_to_k);
    S_from_k = V_from_k .* conj(Ibranch_k);
    S_to_k   = V_to_k   .* conj(-Ibranch_k);
    LineFlows(k, :) = real(S_from_k); % Active power flow from from_bus to to_bus
    BranchLosses(k, :) = real(S_from_k - S_to_k);
end

%==============================Plotting====================================
% Create Figure 1
figure;

% Subplot 1: Generator Frequencies
subplot(3,1,1);
plot(t, w/(2*pi), 'LineWidth', 2);
grid on;
set(gca,'FontSize',12);
legend({'G1','G2','G3'}, 'FontSize',10);
xlabel('Time (sec.)','FontSize',12,'FontWeight','bold');
ylabel('Frequency (Hz)','FontSize',12,'FontWeight','bold');
title('Generator Frequencies Over Time','FontSize',14,'FontWeight','bold');

% Subplot 2: Generator Active Power Outputs
subplot(3,1,2);
plot(t, Pgen, 'LineWidth', 2);
grid on;
set(gca,'FontSize',12);
legend({'G1','G2','G3'}, 'FontSize',10);
xlabel('Time (sec.)','FontSize',12,'FontWeight','bold');
ylabel('Active Power Output (pu)','FontSize',12,'FontWeight','bold');
title('Generator Active Power Outputs Over Time','FontSize',14,'FontWeight','bold');

% Subplot 3: Generator Reactive Power Outputs
subplot(3,1,3);
plot(t, Qgen, 'LineWidth', 2);
grid on;
set(gca,'FontSize',12);
legend({'G1','G2','G3'}, 'FontSize',10);
xlabel('Time (sec.)','FontSize',12,'FontWeight','bold');
ylabel('Reactive Power Output (pu)','FontSize',12,'FontWeight','bold');
title('Generator Reactive Power Outputs Over Time','FontSize',14,'FontWeight','bold');


% Create Figure 2
figure;

% Subplot 1: Bus Voltages
subplot(3,1,1);
plot(t, V, 'LineWidth', 2);
grid on;
set(gca,'FontSize',12);
legend(arrayfun(@(x) sprintf('Bus %d', x), 1:n, 'UniformOutput', false), 'FontSize',8);
xlabel('Time (sec.)','FontSize',12,'FontWeight','bold');
ylabel('Voltage Magnitude (pu)','FontSize',12,'FontWeight','bold');
title('Bus Voltages Over Time','FontSize',14,'FontWeight','bold');

% Subplot 2: Line Flows
subplot(3,1,2);
plot(t, LineFlows, 'LineWidth', 2);
grid on;
set(gca,'FontSize',12);
legend(arrayfun(@(x) sprintf('Line %d', x), 1:num_branches, 'UniformOutput', false), 'FontSize',8);
xlabel('Time (sec.)','FontSize',12,'FontWeight','bold');
ylabel('Active Power Flow (pu)','FontSize',12,'FontWeight','bold');
title('Line Flows Over Time','FontSize',14,'FontWeight','bold');

% Subplot 3: Total Losses
subplot(3,1,3);
plot(t, TotalLosses, 'LineWidth', 2);
grid on;
set(gca,'FontSize',12);
xlabel('Time (sec.)','FontSize',12,'FontWeight','bold');
ylabel('Total Losses (pu)','FontSize',12,'FontWeight','bold');
title('Total Losses Over Time','FontSize',14,'FontWeight','bold');

