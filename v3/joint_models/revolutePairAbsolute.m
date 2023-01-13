classdef revolutePairAbsolute
    %revoluteJointWithRotor 
    
    properties
        nv = 2
        nq = 2
        jointAxis
        XtreeInternal
        output_body = 2
        bodies = 2
    end
    
    methods        
        function [Xup, S, Sd, v,derivs] = kinematics(obj, Xtree, q, qdot, vp)
            % Generalized coordaintes are the absolute angles - o/w same as
            % two revolutes.
            
            
            [XJ{1} , S1 ] = jcalc( ['R' obj.jointAxis{1}], q(1) );
            [XJ{2}, S2 ] = jcalc( ['R' obj.jointAxis{2}], q(2) - q(1) );
            
            X1p = XJ{1}*Xtree(1:6,:);
            X21 = XJ{2}*obj.XtreeInternal;
            
            Xup(1:6,:) = X1p;
            Xup(7:12,:)= X21*X1p;
            
            S  = [S1 zeros(6,1); X21*S1-S2 S2];
            
            if nargout > 2
                v  = Xup*vp + S*qdot;
                v1 = v(1:6);
                v2 = v(7:12);
                
                S1d = crm(v1)*S1;
                S2d = crm(v2)*S2;
                
                Sd = [S1d zeros(6,1) ; X21*S1d-S2d S2d];
            end
            if nargout > 4
                % Sdot (of the open dot variety) lower left corner =
                %                           -crm(S2)*X21*S1*(qdot2 - qdot1)
                % So lower part of Sdot*qdot = 
                %                   -crm(S2)*X21*S1*(qdot2 - qdot1)*qdot1
                
                derivs.Sdotqd_q = zeros(12,2);
                derivs.Sdotqd_qd= zeros(12,2);
                
                dX21_dq1 = crm(S2)*X21;
                dX21_dq2 = -crm(S2)*X21;
                
                derivs.Sdotqd_qd(7:12,1) =  -crm(S2)*X21*S1*(qdot(2) -2*qdot(1));
                derivs.Sdotqd_qd(7:12,2) =  -crm(S2)*X21*S1*(qdot(1)           );
                
                derivs.Sdotqd_q(7:12,1) =  -crm(S2)*dX21_dq1*S1*(qdot(2) - qdot(1))*qdot(1);
                derivs.Sdotqd_q(7:12,2) =  -crm(S2)*dX21_dq2*S1*(qdot(2) - qdot(1))*qdot(1);
            end 
        end
    end
end

