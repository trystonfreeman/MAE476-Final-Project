clear all
clc

syms a mu Omega i u t J2 R
%% Generic
rPQW = [a 0 0]';
vPQW = [0 sqrt(mu/a) 0]';

R1 = [cos(-Omega) sin(-Omega) 0;
      -sin(-Omega) cos(-Omega) 0;
      0 0 1];
R2 = [1 0 0;
    0 cos(-i) sin(-i);
    0 -sin(-i) cos(-i)];
R3 = [cos(-u) sin(-u) 0;
     -sin(-u) cos(-u) 0;
     0 0 1];

R_PQW = R1*R2*R3;

r = R_PQW*rPQW;
v = R_PQW*vPQW;

r = subs(r,[u Omega],[u+2*pi()*t/(2*pi()*sqrt(a^3/mu)), Omega-3/2*J2*sqrt(mu/a^3)*(R/a)^2*cos(i)*t]);
v = subs(v,[u Omega],[u+2*pi()*t/(2*pi()*sqrt(a^3/mu)), Omega-3/2*J2*sqrt(mu/a^3)*(R/a)^2*cos(i)*t]);

h = cross(r,v);

%% Satilite 1
syms Omega1 i1 u1
r1 = subs(r,[u Omega i],[u1 Omega1 i1]);
v1 = subs(v,[u Omega i],[u1 Omega1 i1]);
h1 = subs(h,[u Omega i],[u1 Omega1 i1]);

%% Satilite 2
syms Omega2 i2 u2
r2 = subs(r,[u Omega i],[u2 Omega2 i2]);
v2 = subs(v,[u Omega i],[u2 Omega2 i2]);
h2 = subs(h,[u Omega i],[u2 Omega2 i2]);

%% Intersection Point
N = cross(h1,h2)/norm(cross(h1,h2));

r_int1 = a*N;
r_int2 = -a*N;

distance1 = norm(r1-r_int1);
distance2 = norm(r1-r_int2);

%% Plot
clc
T = 2*pi*sqrt(7378.14^3/398600.44);
t = 0:10:(2*pi*sqrt(7378.14^3/398600.44));

r1t = matlabFunction(subs(r1, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));
v1t = matlabFunction(subs(v1, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));
h1t = matlabFunction(subs(h1, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));

d1t = matlabFunction(subs(distance1, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));

r2t = matlabFunction(subs(r2, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));


pos1 = r1t(t)';
vel1 = v1t(t)';
h1 = h1t(t)';

d1 = d1t(t);

pos2 = r2t(t)';

%%

clc
syms t
sdf1 = (subs(distance1, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));
sdf2 = (subs(distance2, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));

intTime1 = double(vpasolve(sdf1 == 0,t,[0 6400]));
intTime2 = double(vpasolve(sdf2 == 0,t,[0 6400]));
intPos = r1t(min(intTime1,intTime2));

%% Figures
close all
figure(1)
hold on
plot3(pos1(:,1),pos1(:,2),pos1(:,3))
plot3(pos2(:,1),pos2(:,2),pos2(:,3))
plot3(intPos(1),intPos(2),intPos(3),'o')
plot3(pos1(1,1),pos1(1,2),pos1(1,3),'o')
quiver3(pos1(1,1),pos1(1,2),pos1(1,3),1000*vel1(1,1),1000*vel1(1,2),1000*vel1(1,3))

hold off

%% Function ---------------------------------------------------------------
clear
close all
clc

sat1.a = 7378.14;
sat1.i = 80/180*pi;
sat1.Omega = 120/180*pi;
sat1.u = 240/180*pi;

sat2.a = 7378.14;
sat2.i = 80/180*pi;
sat2.Omega = 240/180*pi;
sat2.u = 120/180*pi;

timestart = 0;

[dv,impulsetime] = PlaneChange(sat1,sat2,timestart)

%% First function
function [deltav,t_imp] = PlaneChange(sat1,sat2,timestart)
% Constants
J2 = 0.0010826269;
R = 6378.14;
mu = 398600.44;

% Inputs
a = sat1.a;

Omega1 = sat1.Omega;
u1 = sat1.u;
i1 = sat1.i;

Omega2 = sat2.Omega;
u2 = sat2.u;
i2 = sat2.i;

T = 2*pi*sqrt(a^3/mu);

syms t

r1 = [conj(a)*(cos(u1 + t/(a^3/mu)^(1/2))*cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2)) - sin(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*cos(i1));
    conj(a)*(cos(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2)) + cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*sin(u1 + t/(a^3/mu)^(1/2))*cos(i1));
    sin(u1 + t/(a^3/mu)^(1/2))*conj(a)*sin(i1)];
v1 = [-conj((mu/a)^(1/2))*(cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*sin(u1 + t/(a^3/mu)^(1/2)) + cos(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*cos(i1));
      -conj((mu/a)^(1/2))*(sin(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2)) - cos(u1 + t/(a^3/mu)^(1/2))*cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*cos(i1));
      conj((mu/a)^(1/2))*cos(u1 + t/(a^3/mu)^(1/2))*sin(i1)];
h1 = [conj((mu/a)^(1/2))*cos(u1 + t/(a^3/mu)^(1/2))*conj(a)*sin(i1)*(cos(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2)) + cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*sin(u1 + t/(a^3/mu)^(1/2))*cos(i1)) + conj((mu/a)^(1/2))*sin(u1 + t/(a^3/mu)^(1/2))*conj(a)*sin(i1)*(sin(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2)) - cos(u1 + t/(a^3/mu)^(1/2))*cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*cos(i1));
    - conj((mu/a)^(1/2))*cos(u1 + t/(a^3/mu)^(1/2))*conj(a)*sin(i1)*(cos(u1 + t/(a^3/mu)^(1/2))*cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2)) - sin(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*cos(i1)) - conj((mu/a)^(1/2))*sin(u1 + t/(a^3/mu)^(1/2))*conj(a)*sin(i1)*(cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*sin(u1 + t/(a^3/mu)^(1/2)) + cos(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*cos(i1));
    conj((mu/a)^(1/2))*conj(a)*(cos(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2)) + cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*sin(u1 + t/(a^3/mu)^(1/2))*cos(i1))*(cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*sin(u1 + t/(a^3/mu)^(1/2)) + cos(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*cos(i1)) - conj((mu/a)^(1/2))*conj(a)*(cos(u1 + t/(a^3/mu)^(1/2))*cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2)) - sin(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*cos(i1))*(sin(u1 + t/(a^3/mu)^(1/2))*sin(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2)) - cos(u1 + t/(a^3/mu)^(1/2))*cos(Omega1 - (3*J2*R^2*t*cos(i1)*(mu/a^3)^(1/2))/(2*a^2))*cos(i1))];

r2 = [conj(a)*(cos(u2 + t/(a^3/mu)^(1/2))*cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2)) - sin(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*cos(i2));
    conj(a)*(cos(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2)) + cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*sin(u2 + t/(a^3/mu)^(1/2))*cos(i2));
    sin(u2 + t/(a^3/mu)^(1/2))*conj(a)*sin(i2)];
v2 = [-conj((mu/a)^(1/2))*(cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*sin(u2 + t/(a^3/mu)^(1/2)) + cos(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*cos(i2));
    -conj((mu/a)^(1/2))*(sin(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2)) - cos(u2 + t/(a^3/mu)^(1/2))*cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*cos(i2));
    conj((mu/a)^(1/2))*cos(u2 + t/(a^3/mu)^(1/2))*sin(i2)];
h2 = [conj((mu/a)^(1/2))*cos(u2 + t/(a^3/mu)^(1/2))*conj(a)*sin(i2)*(cos(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2)) + cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*sin(u2 + t/(a^3/mu)^(1/2))*cos(i2)) + conj((mu/a)^(1/2))*sin(u2 + t/(a^3/mu)^(1/2))*conj(a)*sin(i2)*(sin(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2)) - cos(u2 + t/(a^3/mu)^(1/2))*cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*cos(i2));
    - conj((mu/a)^(1/2))*cos(u2 + t/(a^3/mu)^(1/2))*conj(a)*sin(i2)*(cos(u2 + t/(a^3/mu)^(1/2))*cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2)) - sin(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*cos(i2)) - conj((mu/a)^(1/2))*sin(u2 + t/(a^3/mu)^(1/2))*conj(a)*sin(i2)*(cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*sin(u2 + t/(a^3/mu)^(1/2)) + cos(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*cos(i2));
    conj((mu/a)^(1/2))*conj(a)*(cos(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2)) + cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*sin(u2 + t/(a^3/mu)^(1/2))*cos(i2))*(cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*sin(u2 + t/(a^3/mu)^(1/2)) + cos(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*cos(i2)) - conj((mu/a)^(1/2))*conj(a)*(cos(u2 + t/(a^3/mu)^(1/2))*cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2)) - sin(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*cos(i2))*(sin(u2 + t/(a^3/mu)^(1/2))*sin(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2)) - cos(u2 + t/(a^3/mu)^(1/2))*cos(Omega2 - (3*J2*R^2*t*cos(i2)*(mu/a^3)^(1/2))/(2*a^2))*cos(i2))];

N = cross(h1,h2)/norm(cross(h1,h2));

r_int1 = a*N;
r_int2 = -a*N;

distance1 = norm(r1-r_int1);
distance2 = norm(r1-r_int2);

distance1 = matlabFunction(distance1);
distance2 = matlabFunction(distance2);

intTime1 = fsolve(distance1,timestart + T/2);
intTime2 = fsolve(distance2,timestart + T/2);

t_imp = min([intTime1,intTime2]);

v1_imp = double(subs(v1,t_imp));
v2_imp = double(subs(v2,t_imp));

deltav = v2_imp - v1_imp;

end




