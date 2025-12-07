function [dt,dv] = Phase(sat1,sat2)

    phase_angle = sat2.u - sat1.u;
    dt = sat2.T* (1 + phase_angle/360);
    a_phase = (sat1.mu * (dt/2*pi)^2)^(1/3);
    dv = 2* abs(sqrt(sat2.mu*(2/sat2.a - 1/a_phase)) - sqrt(sat2.mu/sat2.a));
end