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
        u_dot
        r
        v
        n
        T
        h
    end
    methods
        function obj = satellite(a,i,omega,u)
            obj.a = a;
            obj.i = i;
            obj.omega = omega;
            obj.u = u;
            
            [obj.r,obj.v] = get_IJK(obj);
            obj.h = cross(obj.r,obj.v);
            obj.n = sqrt(obj.mu/obj.a^3);
            obj.omega_dot = -1.5 * obj.J2*obj.n*(obj.R/obj.a)^2*cosd(obj.i);
            obj.T = 2*pi*sqrt(obj.a^3/obj.mu);
            obj.u_dot =360/obj.T;
        end
        function obj = propagate(obj,dt)

            obj.u = obj.u + dt*obj.u_dot;
            obj.omega = obj.omega + (obj.omega_dot)*(180/pi)*dt;
            [obj.r,obj.v] = get_IJK(obj);
            obj.h = cross(obj.r,obj.v);
        end

        function [r,v] = get_IJK(obj)
            R3_omega = [cosd(-obj.omega), sind(-obj.omega), 0;...
                        sind(-obj.omega),-cosd(-obj.omega), 0;...
                        0              , 0                1];
            R3_u = [cosd(-obj.u), sind(-obj.u), 0;...
                    sind(-obj.u),-cosd(-obj.u), 0;...
                    0          , 0          , 1];
            R1_i = [1 0 0;...
                    0 cosd(-obj.i), sind(-obj.i);...
                    0 sind(-obj.i),-cosd(-obj.i)];
            r = R3_omega*R1_i*R3_u*[obj.a;0;0];
            r = r*obj.a/norm(r); % effectively normalizes transformation
            v = R3_omega*R1_i*R3_u*[0; sqrt(obj.mu/obj.a); 0];
            v = sqrt(obj.mu/obj.a)*v/norm(v);
        end
    end
end