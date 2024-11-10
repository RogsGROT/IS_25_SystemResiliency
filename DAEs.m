function out = DAEs(t,x)
global Vref ws PC Ki Yabs Yang Dm PL PL0 QL m n H Xd Xdp Xdpp Xq Xqp Xqpp Td0p ...
       Td0pp Tq0p Tq0pp Rs Xls KA TA KE TE KF TF Ax Bx TCH TSV RD MH V0 alpha beta bb1 bb2 aa1 aa2

%==================Defining the states and variables=======================
Eqp   = x(1:m);
Si1d  = x(m+1:2*m);
Edp   = x(2*m+1:3*m);
Si2q  = x(3*m+1:4*m);
Delta = x(4*m+1:5*m);
w     = x(5*m+1:6*m);
Efd   = x(6*m+1:7*m);
RF    = x(7*m+1:8*m);
VR    = x(8*m+1:9*m);
TM    = x(9*m+1:10*m);
PSV   = x(10*m+1:11*m);
xi    = x(11*m+1:12*m); % AGC integrator state variable
Id    = x(12*m+1:13*m);
Iq    = x(13*m+1:14*m);
V     = x(14*m+1:14*m+n);
TH    = x(14*m+n+1:14*m+2*n);
COI   = sum(w.*MH')/sum(MH); % Center of Inertia frequency

%===========================Load Increase Over Short Interval=========================

% Increase PL from t = 5s to t = 15s to reach maximum generator output
k = 7; % Scaling factor calculated to match maximum generation capacity

if t >= 5 && t <= 15
    % Linear ramp from PL0 to k*PL0 over 10 seconds
    PLt = PL0 + (t - 5)*(k * PL0 - PL0)/10;
elseif t > 15
    PLt = k * PL0;
else
    PLt = PL0;
end

%===============================DAEs=======================================
VG = V(1:m);
PL2 = PLt'.*((V'./V0).^alpha); % Voltage-dependent active loads
QL2 = QL'.*((V'./V0).^beta);   % Voltage-dependent reactive loads

% Generator dynamic equations
S1 = (1./Td0p).*(-Eqp' - (Xd - Xdp).*(Id' - ((Xdp - Xdpp)./(Xdp - Xls).^2).*(Si1d' + (Xdp - Xls).*Id' - Eqp')) + Efd');
S2 = (1./Td0pp).*(-Si1d' + Eqp' - (Xdp - Xls).*Id');
S3 = (1./Tq0p).*(-Edp' + (Xq - Xqp).*(Iq' - ((Xqp - Xqpp)./(Xqp - Xls).^2).*(Si2q' + (Xqp - Xls).*Iq' + Edp')));
S4 = (1./Tq0pp).*(-Si2q' - Edp' - (Xqp - Xls).*Iq');
S5 = w' - COI;
S6 = (ws./(2*H)).*(TM' - ((Xdpp - Xls)./(Xdp - Xls)).*Eqp'.*Iq' - ((Xdp - Xdpp)./(Xdp - Xls)).*Si1d'.*Iq' - ...
    ((Xqpp - Xls)./(Xqp - Xls)).*Edp'.*Id' + ((Xqp - Xqpp)./(Xqp - Xls)).*Si2q'.*Id' - (Xqpp - Xdpp).*Id'.*Iq' - Dm.*(w' - ws));

% Excitation system equations
S7 = (1./TE).*((-(KE + Ax.*exp(Bx.*Efd'))).*Efd' + VR');
S8 = (1./TF).*(-RF' + (KF./TF).*Efd');
S9 = (1./TA).*(-VR' + KA.*RF' - ((KA.*KF)./TF).*Efd' + KA.*(Vref - VG'));

% Turbine and governor equations with AGC
PCi = PC - Ki.*xi'; % Adjusted power command with AGC
S10 = (1./TCH).*(-TM' + PSV');
S11 = (1./TSV).*(-PSV' + PCi - (1./RD).*(w'./ws - 1));
S12 = w'./ws - 1; % AGC integrator: frequency deviation

% Generator stator and power flow equations
Vectorized_angle1 = (TH(bb1) - TH(aa1) - Yang(1:m,1:n));
Vectorized_mag1 = V(1:m) .* V(1:n)' .* Yabs(1:m,1:n);
sum1 = sum(Vectorized_mag1 .* cos(Vectorized_angle1), 2);
sum2 = sum(Vectorized_mag1 .* sin(Vectorized_angle1), 2);

VG = V(1:m); THG = TH(1:m); Angle_diff = Delta - THG;
SE1 = Rs.*Id' - Xqpp.*Iq' - ((Xqpp - Xls)./(Xqp - Xls)).*Edp' + ((Xqp - Xqpp)./(Xqp - Xls)).*Si2q' + VG'.*sin(Angle_diff)';
SE2 = Rs.*Iq' + Xdpp.*Id' - ((Xdpp - Xls)./(Xdp - Xls)).*Eqp' - ((Xdp - Xdpp)./(Xdp - Xls)).*Si1d' + VG'.*cos(Angle_diff)';

PV1 = (Id'.*VG'.*sin(Angle_diff') + Iq'.*VG'.*cos(Angle_diff')) - PL2(1:m) - sum1';
PV2 = (Id'.*VG'.*cos(Angle_diff') - Iq'.*VG'.*sin(Angle_diff')) - QL2(1:m) - sum2';

% Non-generator power flow equations
Vectorized_angle2 = (TH(bb2) - TH(aa2) - Yang(m+1:n,1:n));
Vectorized_mag2 = V(m+1:n) .* V(1:n)' .* Yabs(m+1:n,1:n);
sum3 = sum(Vectorized_mag2 .* cos(Vectorized_angle2), 2);
sum4 = sum(Vectorized_mag2 .* sin(Vectorized_angle2), 2);

PQ1 = -PL2(m+1:n) - sum3';
PQ2 = -QL2(m+1:n) - sum4';

%--------------------------------------------------------------
out = [S1'; S2'; S3'; S4'; S5'; S6'; S7'; S8'; S9'; S10'; S11'; S12'; SE1'; SE2'; PV1'; PQ1'; PV2'; PQ2'];
end
