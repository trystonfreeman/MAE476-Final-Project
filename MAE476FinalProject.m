clc
clear
close all

mu = 1; % CHANGE ME
a_inner = 7378.14; % [km]
a_outer = 8378.14; % [km]
i_inner = 80; % [deg]
i_outer = 60; % [deg]

inner_sats = satellite.empty(0,6);
inner_sats(1) = satellite(mu,a_inner,i_inner,0,0);
inner_sats(2) = satellite(mu,a_inner,i_inner,0,180);
inner_sats(3) = satellite(mu,a_inner,i_inner,120,60);
inner_sats(4) = satellite(mu,a_inner,i_inner,120,240);
inner_sats(5) = satellite(mu,a_inner,i_inner,240,120);
inner_sats(6) = satellite(mu,a_inner,i_inner,240,300);

outer_sats = satellite.empty(0,3);
outer_sats(1) = satellite(mu,a_outer,i_outer,0,0);
outer_sats(2) = satellite(mu,a_outer,i_outer,120,0);
outer_sats(3) = satellite(mu,a_outer,i_outer,240,0);

servicer_1 = inner_sats(1);
servicer_2 = inner_sats(2);

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

% Maneuver 1

% Maneuver 2

% Maneuver 3

% Maneuver 4

% Maneuver 5

% Maneuver 6

% Maneuver 7

% Maneuver 8