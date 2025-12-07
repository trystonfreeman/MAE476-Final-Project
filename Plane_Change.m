function [deltav,t_imp] = Plane_Change(sat1,sat2,timestart)
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