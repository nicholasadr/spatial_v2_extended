classdef SO3Joint

    properties
        nv = {3}
        nq = {9}
        Xtree
    end
    
    methods        
        function [Xup, S, Sd, v, derivs] = kinematics(obj, q, qdot, vp)
            [XJ , S ] = jcalc( 'SO3', q);            
            Xup = XJ*obj.Xtree;
            if nargout > 2
                v  = Xup*vp + S*qdot;
                Sd = crm(v)*S;
            end
            if nargout > 4
                derivs.Sdotqd_q = zeros(3,1);
                derivs.Sdotqd_qd= zeros(3,1);
            end
        end
    end
end

