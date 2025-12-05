classdef satellite
    properties
        mu
        a
        i
        omega
        u
        r
        r_0
        v_PQW
        v_0_PQW
    end
    methods
        function obj = satellite(mu,a,i,omega,u)
            obj.mu = mu;
            obj.a = a;
            obj.i = i;
            obj.omega = omega;
            obj.u = u;
            
        end
        function propagate(obj)
            
        end
        function maneuver(obj)
            
        end
    end
end