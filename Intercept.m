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
end