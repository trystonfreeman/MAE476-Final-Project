clear
clc

%% Generic 
syms Omega omega mu a i J2 t R u
r = [a*cos(Omega)*cos(omega) - a*sin(Omega)*cos(i)*sin(omega);
a*sin(Omega)*cos(omega) + a*cos(Omega)*cos(i)*sin(omega);
sin(i)*sin(omega)];
% v = [-sqrt(mu/a)*cos(Omega)*sin(omega) - sqrt(mu/a)*sin(Omega)*cos(i)*cos(omega);
% sqrt(mu/a)*cos(Omega)*cos(i)*cos(omega) - sqrt(mu/a)*sin(Omega)*sin(omega);
% sqrt(mu/a)*cos(omega)*sin(i)];

v = [-(cos(Omega)*sin(omega) + sin(Omega)*cos(i)*cos(omega))*(mu/a)^(1/2);
-(sin(Omega)*sin(omega) - cos(Omega)*cos(i)*cos(omega))*(mu/a)^(1/2);
                                      cos(omega)*sin(i)*(mu/a)^(1/2)];

h = cross(r,v);
h_u = subs(h,omega,u);

%h_subs = subs(h_u,[u Omega],[u+360*t/(2*pi()*sqrt(a^3/mu)), Omega-3/2*J2*sqrt(mu/a^3)*(R/a)^2*cos(i)*t]);
h_subs = subs(h_u,[u Omega],[u+2*pi()*t/(2*pi()*sqrt(a^3/mu)), Omega-3/2*J2*sqrt(mu/a^3)*(R/a)^2*cos(i)*t]);

%% for satalite 1
syms Omega1 i1 u1
h1 = subs(h_subs,[Omega i u],[Omega1 i1 u1]);

r1 = subs(r,omega,u);

r1 = subs(r1,[u Omega],[u+2*pi()*t/(2*pi()*sqrt(a^3/mu)), Omega-3/2*J2*sqrt(mu/a^3)*(R/a)^2*cos(i)*t]);

r1 = subs(r1,[Omega i u],[Omega1 i1 u1]);

v1 = subs(v,omega,u);

v1 = subs(v1,[u Omega],[u+2*pi()*t/(2*pi()*sqrt(a^3/mu)), Omega-3/2*J2*sqrt(mu/a^3)*(R/a)^2*cos(i)*t]);

v1 = subs(v1,[Omega i u],[Omega1 i1 u1]);

%% for satalite 2
syms Omega2 i2 u2
h2 = subs(h_subs,[Omega i u],[Omega2 i2 u2]);
r2 = subs(r1,[Omega1 i1 u1],[Omega2 i2 u2]);

v2 = subs(r2,[Omega1 i1 u1],[Omega2 i2 u2]);

%% Calculate intersetion
n = cross(h1,h2)/(norm(cross(h1,h2)));

r_pc(:,1) = a.*n;
r_pc(:,2) = -a.*n;

difffunction1 = (r1-r_pc(:,1));
difffunction2 = (r1-r_pc(:,2));

%% 
interfun1 = norm(subs(difffunction1, ...
      [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], ...
      [2*pi/3 ...
       4*pi/9 ...
       4*pi/3 ...
       4*pi/3 ...
       4*pi/9 ...
       2*pi/3 ...
       7378.14 ...
       398600.44 ...
       0.0010826269 ...
       6378.14]));

interfun2 = norm(subs(difffunction2, ...
      [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], ...
      [2*pi/3 ...
       4*pi/9 ...
       4*pi/3 ...
       4*pi/3 ...
       4*pi/9 ...
       2*pi/3 ...
       7378.14 ...
       398600.44 ...
       0.0010826269 ...
       6378.14]));

interfun1 = matlabFunction(interfun1);
interfun2 = matlabFunction(interfun2);

T = 2*pi*sqrt(7378.14^3/398600.44);

t = 0:10:T;
%% Time functions
r1t = matlabFunction(subs(r1, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));
r2t = matlabFunction(subs(r2, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));

v1t = matlabFunction(subs(v1, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));
v2t = matlabFunction(subs(v2, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));


h1t = matlabFunction(subs(h1, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));
h2t = matlabFunction(subs(h2, [Omega1 i1 u1 Omega2 i2 u2 a mu J2 R], [2*pi/3 4*pi/9 4*pi/3 4*pi/3 4*pi/9 2*pi/3 7378.14 398600.44 0.0010826269 6378.14]));
%% Plots
close all
t = 1:10:6400;
r1vals = r1t(t)';
r2vals = r2t(t)';

v1vals = v1t(t)';
v2vals = v2t(t)';

h1vals = h1t(t)';
h2vals = h2t(t)';

figure(1)
hold on
plot(t,h1vals(:,1),t,h1vals(:,2),t,h1vals(:,3))

figure(2)
plot(t,r1vals(:,1),t,r1vals(:,2),t,r1vals(:,3))

figure(3)
plot(t,v1vals(:,1),t,v1vals(:,2),t,v1vals(:,3))

m = 100;
figure(4)
hold on
plot3(r1vals(:,1),r1vals(:,2),r1vals(:,3))
quiver3(r1vals(1:m,1),r1vals(1:m,2),r1vals(1:m,3),v1vals(1:m,1),v1vals(1:m,2),v1vals(1:m,3))
hold off

%%% Faster version
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

