classdef SE3Joint

    properties
        nv = {6}
        nq = {16}
        Xtree
    end
    
    methods        
        function [Xup, S, Sd, v, derivs] = kinematics(obj, q, qdot, vp)
            [XJ , S ] = jcalc( 'SE3', q);            
            Xup = XJ*obj.Xtree;
            if nargout > 2
                v  = Xup*vp + S*qdot;
                Sd = crm(v)*S;
            end
            if nargout > 4
                derivs.Sdotqd_q = zeros(6,1);
                derivs.Sdotqd_qd= zeros(6,1);
            end
        end
    end
end

