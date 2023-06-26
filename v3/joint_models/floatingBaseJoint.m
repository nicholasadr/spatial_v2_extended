classdef floatingBaseJoint

    properties
        nv = 6
        nq = 7
        bodies = 1
        output_body=1
    end
    
    methods        
        function [Xup, S, Sd, v, derivs] = kinematics(obj, Xtree, q, qdot, vp)
            [XJ , S ] = jcalc( 'Fb', q);              
            Xup = XJ*Xtree;
            if nargout > 2
                v  = Xup*vp + S*qdot;
                Sd = crm(v)*S;
            end
            if nargout > 4
                derivs.S_q = zeros(6,1,1);
                derivs.Sdotqd_q = zeros(6,1);
                derivs.Sdotqd_qd= zeros(6,1);
            end
        end
    end
end

