classdef satellite
    properties
        a
        i
        omega
        u
        r
        r_0
        v
        v_0
    end
    methods
        function obj = satellite(a,i,omega,u)
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