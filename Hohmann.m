% Hohmann transfer between 2 circular orbits around a given body
function [dv_1,dv_2,dt,dtheta] = Hohmann(sat1,sat2)
    v_1 = sqrt(sat1.mu/sat1.a);
    v_2 = sqrt(sat2.mu/sat2.a);
    a_t = (sat1.a + sat2.a)/2;

    dv_1 = sqrt(sat1.mu*(2/sat1.a -1/a_t)) - v_1; % initial delta v required (scalar)
    dv_2 = v_2 - sqrt(sat2.mu*(2/sat2.a -1/a_t)); % secondary delta v required (scalar)
    dt = pi*sqrt(a_t^3/sat2.mu); % TOF for transfer (scalar)
    dtheta = 2*pi*(dt/(2*pi*sqrt(sat2.a^3/sat1.mu))); % angle sat 2 sweeps during transfer
end