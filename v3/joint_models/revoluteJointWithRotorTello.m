classdef revoluteJointWithRotorTello
    %revoluteJointWithRotorTello based on revoluteJointWithRotor
    
    properties
        nv = 1
        nq = 1
        rotorAxis
        jointAxis
        gearRatio
        output_body = 2
        bodies = 2
    end
    
    methods        
        function [Xup, S, Sd, v, derivs] = kinematics(obj, Xtree, y, y_dot, vp)
            [XJ_rotor, S_rotor ] = jcalc( obj.rotorAxis, y*obj.gearRatio);
            [XJ_link , S_link ] = jcalc( obj.jointAxis, y);

            Xr = XJ_rotor * Xtree(1:6,:);
            Xl = XJ_link * Xtree(7:12,:);
            Xup = [Xr ; Xl];

            S = [S_rotor*obj.gearRatio ; S_link];

            if nargout > 2
                v  = Xup*vp + S*y_dot;
                Sd = crm(v)*S;
            end
        end
    end
end

