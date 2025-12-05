classdef satellite
    properties
        J2 = 0.0010826269;
        R = 6378.14;
        mu = 398600.44;
        a
        i
        omega
        omega_dot
        u
        r_0_PQW
        r_PQW
        v_PQW
        v_0_PQW
        n
        T
    end
    methods
        function obj = satellite(a,i,omega,u)
            obj.a = a;
            obj.i = i;
            obj.omega = omega;
            obj.u = u;
            
            r_vec = get_PQW(obj);
            obj.r_0_PQW = r_vec(1:3);
            obj.v_0_PQW = r_vec(4:6);
            obj.r_PQW = obj.r_0_PQW;
            obj.v_PQW = obj.v_0_PQW;
            obj.n = sqrt(obj.mu/obj.a^3);
            obj.omega_dot = -1.5 * obj.J2*obj.n*(obj.R/obj.a)^2*cosd(obj.i);
            obj.T = 2*pi*sqrt(obj.a^3/obj.mu);
            
        end
        function obj = propagate(obj,dt)
            obj.u = obj.u+dt/obj.T;
            obj.omega = obj.omega + obj.omega_dot*dt;
        end

        function r_vec = get_PQW(obj)
            v_PQW = sqrt(obj.mu/obj.a) * [-sind(obj.u); cosd(obj.u); 0];
            r_PQW = [obj.a*cosd(obj.u);obj.a*sind(obj.u);0];
            r_vec = [r_PQW;v_PQW];
        end
    end
end