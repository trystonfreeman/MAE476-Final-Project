clc
clear
close all
tic
dt = 1;  % [s]
T_total = 40000; %[s]
t = 0:dt:T_total;
N = T_total/dt + 1;
N_man = 12; % Number of maneuvers

a_inner = 7378.14; % [km]
a_outer = 8378.14; % [km]
i_inner = 80; % [deg]
i_outer = 60; % [deg]
inner_sats = satellite.empty(0,6);
inner_sats(1) = satellite(a_inner,i_inner,0,0);
inner_sats(2) = satellite(a_inner,i_inner,0,180);
inner_sats(3) = satellite(a_inner,i_inner,120,60);
inner_sats(4) = satellite(a_inner,i_inner,120,240);
inner_sats(5) = satellite(a_inner,i_inner,240,120);
inner_sats(6) = satellite(a_inner,i_inner,240,300);


outer_sats = satellite.empty(0,3);
outer_sats(1) = satellite(a_outer,i_outer,0,0);
outer_sats(2) = satellite(a_outer,i_outer,120,0);
outer_sats(3) = satellite(a_outer,i_outer,240,0);

servicer_1 = inner_sats(1);
servicer_2 = inner_sats(1);

inner_sat_pos = NaN(18,N);
outer_sat_pos = NaN(9,N);


%% Simulate Satellites (Assumptions)
for i=1:N

    for j=1:6 % Inner Sats
        % Logs values at each time step
        inner_sat_pos(3*j-2:3*j,i) = inner_sats(j).r;
        
        % propagate next time step
        inner_sats(j) = inner_sats(j).propagate(t(i));
        
    end

    for j = 1:3 % Outer Sats
        outer_sat_pos(3*j-2:3*j,i) = outer_sats(j).r;
        outer_sats(j) = outer_sats(j).propagate(dt);
    end
    
end

%% Simulate Satellites (ODE45)


%% Plan Maneuvers (Option 1)

% THIS SECTION SHOULD BE MOSTLY GOOD
% Servicer 1 (Staying in center)
tserv1_0 = 0;
dv_serv1 = 0;
    % Maneuver 1 (sat 1 -> sat 2)
    % Phase
    [dt1,dv1a,dv1b] = Phase(servicer_1,inner_sats(2));
    tserv1_1 = tserv1_0 + dt1;
    dv_serv1 = dv_serv1 + norm(dv1a) + norm(dv1b);

    servicer_1.phase(inner_sats(2));

    % Maneuver 2-3 (sat 2 -> sat 3)
    % Plane Change
    [tserv1_2,dv2] = Intercept(inner_sats(2),inner_sats(3),tserv1_1);
    dv_serv1 = dv_serv1 + norm(dv2);
    inner_sats(2).propagate(tserv1_2);
    inner_sats(3).propagate(tserv1_2);

    servicer_1.propagate(tserv1_2);
    servicer_1.plane_change(inner_sats(3));

    % Phase
    [dt3,dv3a,dv3b] = Phase(inner_sats(2),inner_sats(3));
    tserv1_3 = tserv1_2 + dt3;
    dv_serv1 = dv_serv1 + norm(dv3a) + norm(dv3b);
    
    % Maneuver 4 (sat 3 -> sat 4)
    % Phase
    [dt4,dv4a,dv4b] = Phase(inner_sats(3),inner_sats(4));
    tserv1_4 = tserv1_3 + dt4;
    dv_serv1 = dv_serv1 + norm(dv4a) + norm(dv4b);

    % Maneuver 5-6 (sat 4 -> sat 5)
    % Plane Change
    [tserv1_5,dv5] = Intercept(inner_sats(4),inner_sats(5),tserv1_4);
    dv_serv1 = dv_serv1 + norm(dv5);
    inner_sats(4).propagate(tserv1_5);
    inner_sats(5).propagate(tserv1_5);
    % Phase
    [dt6,dv6a,dv6b] = Phase(inner_sats(3),inner_sats(4));
    tserv1_6 = tserv1_5 + dt6;
    dv_serv1 = dv_serv1 + norm(dv6a) + norm(dv6b);

    % Maneuver 7 (sat 5 -> sat 6)
    % Phase
    [dt7,dv7a,dv7b] = Phase(inner_sats(3),inner_sats(4));
    tserv1_7 = tserv1_6 + dt7;
    dv_serv1 = dv_serv1 + norm(dv7a) + norm(dv7b);

% THIS NEEDS WORK
% Servicer 2 (Going to outer ring)
    % Maneuver 8-10 (sat 1 -> sat 7)
    % Hohmann
    [dv8a,dv8b,dt8,dtheta8] = Hohmann(servicer_2,outer_sats(1));
    tserv2_1 = dt8;
    dv_serv2 = norm(dv8a) + norm(dv8b);

    % Plane Change
    [tserv2_2,dv9] = Intercept(servicer_2,outer_sats(1),tserv2_1);
    dv_serv2 = dv_serv2 + norm(dv9);

    % Phase
    [dt10,dv10a,dv10b] = Phase(servicer_2,outer_sats(1));
    dv_serv2 = dv_serv2 + norm(dv10a) + norm(dv10b);
    tserv2_3 = tserv2_2 + dt10;

    % Maneuver 11 (sat 7 -> sat 8)
    % Plane Change
    [tserv2_4,dv11] = Intercept(servicer_2,outer_sats(2),tserv2_3);
    dv_serv2 = dv_serv2 + norm(dv11);

    % Maneuver 12 (sat 8 -> sat 9)
    % Plane Change
    [tserv2_5,dv12] = Intercept(servicer_2,outer_sats(3),tserv2_4);
    dv_serv2 = dv_serv2 + norm(dv12);
    dv = dv_serv1 + dv_serv2;

%% Plan Maneuvers (Option 2)

%% Calculate Masses
g0 = 9.81; % [m/s^2]
Isp = 250; %[s]
ms = 300; %[kg]
m_pay = 5; %[kg]
m_final = ms + m_pay;
ve = Isp * g0;

% Maneuver 12
mi = m_final / exp(-dv12*1000/ve);
m_prop12 = mi - m_final;

% Maneuver 11
m_final = m_final + m_prop12;
mi = m_final / exp(-dv11*1000/ve);
m_prop11 = mi - m_final;

% Maneuver 10
m_final = m_final + m_prop11;
mi = m_final / exp(-dv11*1000/ve);
m_prop11 = mi - m_final;


%% Simulate Maneuvers
t_man = NaN(N_man+N,1);
j = 0;
for i = 1:N
    % Servicer 1 Maneuvers
    % Maneuver 1
    t_man(i+j) = t(i); % Overwrites this if a maneuver happens between time steps
    if (i == N)

    elseif (t(i) == 0)
        t_man(i+j) = 0;
        servicer_1.v = servicer_1.v + dv1a;

    elseif ((t(i+1) > tserv1_1) & (t(i) < tserv1_1))
        t_man(i+j) = tserv1_1;
        j = j+1;
        servicer_1.v = servicer_1.v + dv1b;
    
    % Maneuver 2
    elseif ((t(i+1) > tserv1_2) & (t(i) < tserv1_2))
        t_man(i+j) = tserv1_2;
        j = j+1;
        servicer_1.v = servicer_1.v + dv2;
    
    % Maneuver 3
        servicer_1.v = servicer_1.v + dv3a;
    elseif ((t(i+1) > tserv1_3) & (t(i) < tserv1_3))
        t_man(i+j) = tserv1_3;
        j = j+1;
        servicer_1.v = servicer_1.v + dv3b;
    % Maneuver 4
        servicer_1.v = servicer_1.v + dv4a;
    elseif ((t(i+1) > tserv1_4) & (t(i) < tserv1_4))
        t_man(i+j) = tserv1_4;
        j = j+1;
        servicer_1.v = servicer_1.v + dv4b;
    % Maneuver 5
    elseif ((t(i+1) > tserv1_5) & (t(i) < tserv1_5))
        t_man(i+j) = tserv1_5;
        j = j+1;
        servicer_1.v = servicer_1.v + dv5;
    % Maneuver 6
    servicer_1.v = servicer_1.v + dv6a;
    elseif ((t(i+1) > tserv1_6) & (t(i) < tserv1_6))
        t_man(i+j) = tserv1_6;
        j = j+1;
        servicer_1.v = servicer_1.v + dv6b;
    % Maneuver 7
    servicer_1.v = servicer_1.v + dv7a;
    elseif ((t(i+1) > tserv1_7) & (t(i) < tserv1_7))
        t_man(i+j) = tserv1_7;
        j = j+1;
        servicer_1.v = servicer_1.v + dv7b;
    end
    
    % Servicer 2 Maneuvers
    % Maneuver 8
    if (i == N)

    elseif (t(i) == 0)
        t_man(i+j) = 0;
        servicer_2.v = servicer_2.v + dv8a;
    elseif ((t(i+1) > tserv2_1) & (t(i) < tserv2_1))
        t_man(i+j) = tserv2_1;
        j = j+1;
        servicer_2.v = servicer_2.v + dv8b;
    
    % Maneuver 9
    elseif ((t(i+1) > tserv2_2) & (t(i) < tserv2_2))
        t_man(i+j) = tserv2_2;
        j = j+1;
        servicer_2.v = servicer_2.v + dv9;
    % Maneuver 10
    servicer_2.v = servicer_2.v + dv10a;
    elseif ((t(i+1) > tserv2_3) & (t(i) < tserv2_3))
        t_man(i+j) = tserv2_3;
        j = j+1;
        servicer_2.v = servicer_2.v + dv10b;
    % Maneuver 11
    elseif ((t(i+1) > tserv2_4) & (t(i) < tserv2_4))
        t_man(i+j) = tserv2_4;
        j = j+1;
        servicer_2.v = servicer_2.v + dv11;
    % Maneuver 12
    elseif ((t(i+1) > tserv2_5) & (t(i) < tserv2_5))
        t_man(i+j) = tserv2_5;
        j = j+1;
        servicer_2.v = servicer_2.v + dv12;
    end
    
    % Simulate orbit with updated velocity
end
toc

%% Figures
figure("Name","Constellation")
hold on
for i = 1:6
plot3(inner_sat_pos(3*i-2,:),inner_sat_pos(3*i-1,:),inner_sat_pos(3*i,:))
%plot3([0,inner_sat_h(3*i-2,1)],[0,inner_sat_h(3*i-1,1)],[0,inner_sat_h(3*i,1)])
end

for i = 1:3
plot3(outer_sat_pos(3*i-2,:),outer_sat_pos(3*i-1,:),outer_sat_pos(3*i,:))
end
view(3)
figure("Name","Servicer Path")

figure("Name","Pos Error")

%% Functions


% Hohmann transfer between 2 circular orbits around a given body
function [dv_1,dv_2,dt,dtheta] = Hohmann(sat1,sat2)
    v_1 = sqrt(sat1.mu/sat1.a);
    v_2 = sqrt(sat2.mu/sat2.a);
    a_t = (sat1.a + sat2.a)/2;

    dv_1 = sqrt(sat1.mu*(2/sat1.a -1/a_t)) - v_1; % initial delta v required (scalar)
    dv_2 = v_2 - sqrt(sat2.mu*(2/sat2.a -1/a_t)); % secondary delta v required (scalar)
    dt = pi*sqrt(a_t^3/sat2.mu); % TOF for transfer (scalar)
    dtheta = 2*pi*(dt/(2*pi*sqrt(sat2.a^3/sat1.mu))); % angle sat 2 sweeps during transfer
    if (v_1 > v_2)
        direction =  sat1.v/norm(sat1.v);
    else
        direction = -sat1.v/norm(sat1.v);
    end
    dv_1 = dv_1 *   direction; % initial delta v required (vector)
    dv_2 = dv_2 * (-direction); % secondary delta v required (vector)
end



% Intersection of two co-radial circular orbits
function [t_intercept,dv] = Intercept(sat1,sat2,t_0)
    u1 = @(t) wrapTo360(sat1.u_0 + sat1.u_dot*t);
    u2 = @(t) wrapTo360(sat2.u_0 + sat2.u_dot*t);
    omega1 = @(t) sat1.omega_0 + sat1.omega_dot*t;
    omega2 = @(t) sat2.omega_0 + sat2.omega_dot*t;
    R_i = [1 0 0;
           0 cosd(-sat1.i)  sind(-sat1.i);
           0 sind(-sat1.i) -cosd(-sat1.i)];
    R_omega1 = @(t)[cosd(-omega1(t)), sind(-omega1(t)), 0;
                    sind(-omega1(t)),-cosd(-omega1(t)), 0;
                    0 0 1];
    R_omega2 = @(t)[cosd(-omega2(t)), sind(-omega2(t)), 0;
                    sind(-omega2(t)),-cosd(-omega2(t)), 0;
                    0 0 1];
    R_u1 = @(t) [cosd(-u1(t)), sind(-u1(t)), 0;...
                 sind(-u1(t)),-cosd(-u1(t)), 0;...
                 0           , 0           , 1];
    R_u2 = @(t) [cosd(-u2(t)), sind(-u2(t)), 0;...
                 sind(-u2(t)),-cosd(-u2(t)), 0;...
                 0           , 0           , 1];

    r1 = @(t)R_omega1(t)*R_i*R_u1(t)*[sat1.a;0;0];
    r2 = @(t)R_omega2(t)*R_i*R_u2(t)*[sat2.a;0;0];

    v1 = @(t)R_omega1(t)*R_i*R_u1(t)*sqrt(sat1.mu/sat1.a)*[0;1;0];
    v2 = @(t)R_omega2(t)*R_i*R_u2(t)*sqrt(sat2.mu/sat2.a)*[0;1;0];

    h1 = @(t) cross(r1(t),v1(t));
    h2 = @(t) cross(r2(t),v2(t));
    r_intercept = @(t) sat1.a*cross(h1(t),h2(t))/norm(cross(h1(t),h2(t)));
    separation = @(t) min(norm(r_intercept(t) - r1(t)),norm(-r_intercept(t) - r1(t)));

    t_intercept = fsolve(separation,t_0);
    if (t_intercept < t_0)
        t_0 = t_intercept + 0.5*sat1.T;
        t_intercept = fsolve(separation,t_0);
    end
    delta = acosd(cos(sat1.i)^2 +sind(sat1.i)^2*cosd(sat2.omega - sat1.omega));
    dv = 2*sqrt(sat1.mu/sat1.a)*sin(delta/2);

    direction = sign(omega2(t_intercept) - omega1(t_intercept)) * h1(t_intercept)/norm(h1(t_intercept));

    dv = dv * direction;
end



% Phase maneuver between circular orbits
function [dt,dv1,dv2] = Phase(sat1,sat2)

    phase_angle = sat2.u - sat1.u;
    dt = sat2.T* (1 + phase_angle/360);
    a_phase = (sat1.mu * (dt/2*pi)^2)^(1/3);
    dv = 2* abs(sqrt(sat2.mu*(2/sat2.a - 1/a_phase)) - sqrt(sat2.mu/sat2.a));

    direction = sat1.v/norm(sat1.v);

    dv1 = dv * (-direction);
    dv2 = dv * ( direction);
end

% Calculate Mass 
function m_prop = calmass(m_i, m_final, dv)

end