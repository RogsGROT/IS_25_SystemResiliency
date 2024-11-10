function mpc = IEEE9Bus
%IEEE9BUS  Power flow data for 9 bus, 3 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Data from:
%       G. W. Stagg and A. H. El-Abiad, "Computer Methods in Power System
%       Analysis", McGraw-Hill, 1968.

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%   bus_i  type    Pd      Qd      Gs      Bs      area    Vm      Va      baseKV  zone    Vmax    Vmin
mpc.bus = [
    1       3       0       0       0       0       1       1.04    0       16.5    1       1.1     0.9;
    2       2       0       0       0       0       1       1.025   0       18.0    1       1.1     0.9;
    3       2       0       0       0       0       1       1.025   0       13.8    1       1.1     0.9;
    4       1       50      30      0       0       1       1       0       230     1       1.1     0.9;
    5       1       20      10      0       0       1       1       0       230     1       1.1     0.9;
    6       1       20      10      0       0       1       1       0       230     1       1.1     0.9;
    7       1       0       0       0       0       1       1       0       230     1       1.1     0.9;
    8       1       0       0       0       0       1       1       0       230     1       1.1     0.9;
    9       1       0       0       0       0       1       1       0       230     1       1.1     0.9;
];

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    1       71.64   27.05   300     -300    1.04    100     1       250     10;
    2       163     6.7     300     -300    1.025   100     1       300     10;
    3       85      -10.9   300     -300    1.025   100     1       270     10;
];

%% branch data
%   fbus    tbus    r           x           b           rateA   rateB   rateC   ratio   angle   status  angmin  angmax
mpc.branch = [
    4       5       0           0.0576      0           250     250     250     0       0       1       -360    360;
    4       6       0           0.092       0           250     250     250     0       0       1       -360    360;
    5       7       0.017       0.092       0           150     150     150     0       0       1       -360    360;
    6       9       0.039       0.17        0           150     150     150     0       0       1       -360    360;
    7       8       0.0085      0.072       0           150     150     150     0       0       1       -360    360;
    8       9       0.0119      0.1008      0           150     150     150     0       0       1       -360    360;
    2       7       0           0.0625      0           250     250     250     0       0       1       -360    360;
    3       9       0           0.0586      0           300     300     300     0       0       1       -360    360;
    1       4       0           0.0576      0           250     250     250     0       0       1       -360    360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%   1   startup  shutdown    n   x1  y1  ... xn  yn
%   2   startup  shutdown    n   c(n-1)  ... c0
mpc.gencost = [
    2   0       0   3   0.02    2   0;
    2   0       0   3   0.0175  1.75    0;
    2   0       0   3   0.0625  1      0;
];

