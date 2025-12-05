clear
close all
clc

% Convert Orbital Elements to XYZ
[s1r, s1v] = OE2rv(7378.14,80,0,0);
[s2r, s2v] = OE2rv(7378.14,80,0,180);
[s3r, s3v] = OE2rv(7378.14,80,120,60);
[s4r, s4v] = OE2rv(7378.14,80,120,240);
[s5r, s5v] = OE2rv(7378.14,80,240,120)
[s6r, s6v] = OE2rv(7378.14,80,240,300);
[s7r, s7v] = OE2rv(8378.14,60,0,0);
[s8r, s8v] = OE2rv(8378.14,60,120,0);
[s9r, s9v] = OE2rv(8378.14,60,240,0);

% Plot Origional Positions
posmat = [s1r s2r s3r s4r s5r s6r s7r s8r s9r];
velmat = [s1v s2v s3v s4v s5v s6v s7v s8v s9v];

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

T1 = 2*pi()*sqrt(7378.14^3/398600.44);
T2 = 2*pi()*sqrt(8378.14^3/398600.44);
orb = 1;

% [t,q1] = ode45(@dqdt3d,0:100:orb*T1,[s1r; s1v],options);
% [t,q2] = ode45(@dqdt3d,0:100:orb*T1,[s2r; s2v],options);
% [t,q3] = ode45(@dqdt3d,0:100:orb*T1,[s3r; s3v],options);
[t,q4] = ode45(@dqdt3d,0:100:orb*T1,[s4r; s4v],options);
[t,q5] = ode45(@dqdt3d,0:100:orb*T1,[s5r; s5v],options);
% [t,q6] = ode45(@dqdt3d,0:100:orb*T2,[s6r; s6v],options);
% [t,q7] = ode45(@dqdt3d,0:100:orb*T2,[s7r; s7v],options);
% [t,q8] = ode45(@dqdt3d,0:100:orb*T2,[s8r; s8v],options);

%openfig("constellation.fig")

Dist = zeros(length(q4),1);
for i = 1:length(q4(:,1))
    D1 = 10000000*ones(length(q5),1);
    for j = 1:length(q5(:,1))
        D1(j) = norm(q4(i,1:3) - q5(j,1:3));
    end
    Dist(i) = min(D1);
end
midpoint = floor(length(Dist)/2);
minDist1 = min(Dist(1:midpoint));
minInd1 = find(Dist == minDist1)
minDist2 = min(Dist(midpoint+1:end));
minInd2 = find(Dist == minDist2)

hold on
plot3([q4(minInd1,1) q4(minInd2,1)],[q4(minInd1,2) q4(minInd2,2)],[q4(minInd1,3) q4(minInd2,3)],'o','Color','k')
plot3(posmat(1,:),posmat(2,:),posmat(3,:),'o','Color','g','LineWidth',2)
quiver3(0,0,0,0,5000,0,"Color",'k')
quiver3(0,0,0,5000,0,0,"Color",'k')
quiver3(0,0,0,0,0,5000,"Color",'k')
% plot3(q1(:,1),q1(:,2),q1(:,3))
% plot3(q2(:,1),q2(:,2),q2(:,3))
% plot3(q3(:,1),q3(:,2),q3(:,3))
plot3(q4(:,1),q4(:,2),q4(:,3))
plot3(q5(:,1),q5(:,2),q5(:,3))
% plot3(q6(:,1),q6(:,2),q6(:,3))
% plot3(q7(:,1),q7(:,2),q7(:,3))
% plot3(q8(:,1),q8(:,2),q8(:,3))
quiver3(posmat(1,:),posmat(2,:),posmat(3,:),velmat(1,:),velmat(2,:),velmat(3,:))
grid on
axis equal
view(130,30)
hold off

function [r, v] = OE2rv(a,i,Omega,u)
mu = 398600.44;
rPQW = [a*cosd(0) a*sind(0) 0]';
vPQW = sqrt(mu/a)*[-sin(0) cos(0) 0]';

R1 = [cosd(-Omega)  sind(-Omega) 0;
      -sind(-Omega) cosd(-Omega) 0;
      0             0            1];
R2 = [1    0      0;
      0 cosd(-i) sind(-i);
      0 -sind(-i) cosd(-i)];
R3 = [cosd(-u) sind(-u) 0;
     -sind(-u) cosd(-u) 0;
            0        0  1];
PQW_IJK = R1*R2*R3;

r = PQW_IJK*rPQW;
v = PQW_IJK*vPQW;

end

function dqdt = dqdt3d(t,q)
r = sqrt(q(1)^2 + q(2)^2 + q(3)^2);
dqdt(1) = q(4);
dqdt(2) = q(5);
dqdt(3) = q(6);

c = 3/2*.0010826269*398600.44*6378.14^2/r^4;

dqdt(4) = -398600.44*q(1)/r^3 + c*q(1)/r*(5*q(3)^2/r^2 - 1);
dqdt(5) = -398600.44*q(2)/r^3 + c*q(2)/r*(5*q(3)^2/r^2 - 1);
dqdt(6) = -398600.44*q(3)/r^3 + c*q(3)/r*(5*q(3)^2/r^2 - 3);

dqdt = dqdt';
end