clc
clear
close all
dt = 10;
T_total = 3600;
t = 0:dt:T_total;
N = T_total/dt + 1;

mu = 1; % CHANGE ME
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

inner_sats_const_omega = inner_sats;
for i = 1:6
inner_sats_const_omega(i).omega_dot = 0;
end
outer_sats = satellite.empty(0,3);
outer_sats(1) = satellite(a_outer,i_outer,0,0);
outer_sats(2) = satellite(a_outer,i_outer,120,0);
outer_sats(3) = satellite(a_outer,i_outer,240,0);

servicer_1 = inner_sats(1);
servicer_2 = inner_sats(1);

inner_sat_pos = NaN(18,N);
outer_sat_pos = NaN(9,N);

inner_sat_omega = NaN(6,N);
inner_sat_h = NaN(18,N);
inner_sat_u = NaN(6,N);
%% Simulate Satellites (Assumptions)
for i=1:N

    for j=1:6
        % Logs values at each time step
        inner_sat_pos(3*j-2:3*j,i) = inner_sats(j).r;
        inner_sat_h(3*j-2:3*j,i) = inner_sats(j).h;
        inner_sat_u(j,i) = inner_sats(j).u;
        inner_sat_omega(j,i) = inner_sats(j).omega;
        
        % propagate next time step
        inner_sats(j) = inner_sats(j).propagate(dt);
    end
    for j = 1:3
        outer_sat_pos(3*j-2:3*j,i) = outer_sats(j).r;
        outer_sats(j) = outer_sats(j).propagate(dt);
    end
end

figure()
hold on
for i = 1:6
plot3(inner_sat_pos(3*i-2,:),inner_sat_pos(3*i-1,:),inner_sat_pos(3*i,:))
end
sat13_intercept = a_inner*normalize(cross(inner_sat_h(1:3,:),inner_sat_h(7:9,:)));
plot3(sat13_intercept(1,:),sat13_intercept(2,:),sat13_intercept(3,:))
legend()
for i = 1:3
%plot3(outer_sat_pos(3*i-2,:),outer_sat_pos(3*i-1,:),outer_sat_pos(3*i,:))
end
view(3)

figure()
plot(t,inner_sat_u(1,:))
figure()


%% Simulate Satellites (ODE45)


%% Plan Maneuvers

% Maneuver 1

% Maneuver 2

% Maneuver 3

% Maneuver 4

% Maneuver 5

% Maneuver 6

% Maneuver 7

% Maneuver 8

%% Simulate Maneuvers
for i = 1:N
% Maneuver 1

% Maneuver 2

% Maneuver 3

% Maneuver 4

% Maneuver 5

% Maneuver 6

% Maneuver 7

% Maneuver 8
end