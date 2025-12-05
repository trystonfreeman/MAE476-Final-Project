clear all
clc

%% Constants
J2 = 0.0010826269;
mu = 398600.44;
R = 6378.14;

%% Givens
a = 7378.14;
i = 80;
Omega0 = 120;
u0 = 240;

%% Calculated Constants
n = sqrt(mu/a^3); % Mean Motion?
Omegadot = -3/2*J2*n*(R/a)^2*cosd(i);
T = 2*pi*sqrt(a^3/mu);
P = a;


%% Equations
Omega = @(t) Omega0 + Omegadot*t;

u = @(t) wrapTo360(u0 + 360*t/T);

%% Time 
t = 100;

%% Velocity and Position

R1 = [cosd(-Omega(t))  sind(-Omega(t)) 0;
      -sind(-Omega(t)) cosd(-Omega(t)) 0;
      0             0            1];
R2 = [1    0      0;
      0 cosd(-i) sind(-i);
      0 -sind(-i) cosd(-i)];
R3 = [cosd(-u(t)) sind(-u(t)) 0;
     -sind(-u(t)) cosd(-u(t)) 0;
            0        0  1];
PQW_IJK = R1*R2*R3;

rPQW = [a 0 0]';
vPQW = sqrt(mu/P)*[0 1 0]';

r = PQW_IJK*rPQW;
v = PQW_IJK*vPQW;


