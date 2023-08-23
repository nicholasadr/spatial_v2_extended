classdef telloDifferential
    % Tello Differential

    properties
        type
        nv = 2
        nq = 2
        jointAxis
        rotorAxis

        gearRatio = 6

        XtreeInternal
        output_body = 4
        bodies = 4
    end

    methods
        function [Xup, S, Sd, v] = kinematics(obj, Xtree, ql, ql_dot, vp)
            
            % Inverse kinematics to obtain
            % y: post-gearbox rotor angle (minimal coordinate)
            y_fn = casadi.external(strcat(obj.type,'_IK_pos'),...
                                    strcat('./models/casadi_fn/gen_',...
                                           obj.type,...
                                           '_IK_ql_to_y.so'));
            y = full( y_fn(ql) );

            % qr: pre-gearbox rotor angle (independent maximal coordinates)
            qr = obj.gearRatio * y;
            
            [XR1 , Sr1 ] = jcalc( obj.rotorAxis{1}, qr(1) );
            [XR2 , Sr2 ] = jcalc( obj.rotorAxis{2}, qr(2) );
            [XJ1 , S1 ] = jcalc( obj.jointAxis{1}, ql(1) );
            [XJ2 , S2 ] = jcalc( obj.jointAxis{2}, ql(2) );

            Xr1p = XR1 * Xtree(1:6,:);
            Xr2p = XR2 * Xtree(7:12,:);
            X1p = XJ1 * Xtree(13:18,:);
            X21 = XJ2 * obj.XtreeInternal;
            X2p = X21 * X1p;

            Xup = [Xr1p ; Xr2p ; X1p ; X2p];

            G_fn = casadi.external(strcat(obj.type,'_G'),...
                   strcat('./models/casadi_fn/gen_',...
                          obj.type,'_G.so'));
            G = full( G_fn(qr,ql) );

            z = zeros(6,1);
            S = [obj.gearRatio*Sr1              z;
                 z                          obj.gearRatio*Sr2;
                 G(3,1)*S1                  G(3,2)*S1;
                 G(3,1)*X21*S1+G(4,1)*S2    G(3,2)*X21*S1+G(4,2)*S2
                ];

            if nargout > 2
                y_dot_fn = casadi.external(...
                            strcat(obj.type,'_IK_vel'),...
                            strcat('./models/casadi_fn/gen_',...
                                   obj.type,'_IK_dql_to_dy.so'));
                y_dot = full( y_dot_fn(ql,ql_dot) );

                qr_dot = obj.gearRatio * y_dot;

                G_dot_fn = casadi.external(...
                           strcat(obj.type,'_G_dot'),...
                           strcat('./models/casadi_fn/gen_',...
                           obj.type,'_G_dot.so'));
                G_dot = full( G_dot_fn(qr,ql,qr_dot,ql_dot) );

                v = Xup*vp + S*y_dot;
                vr1 = v(1:6);
                vr2 = v(7:12);
                v1 = v(13:18);
                v2 = v(19:24);

                X21S1 = X21*S1;
                v2_relx = crm(S2*ql_dot(2));
                Sr1x = crm(vr1)*Sr1;
                Sr2x = crm(vr2)*Sr2;
                S1x = crm(v1)*S1;
                S2x = crm(v2)*S2;
                X21S1x = crm(v2)*X21S1;

                S_ring = [z                        z;
                          z                        z;
                          G_dot(3,1)*S1            G_dot(3,2)*S1;
                          G_dot(3,1)*X21S1+...
                          G_dot(4,1)*S2+...
                          G(3,1)*(-v2_relx*X21S1)  G_dot(3,2)*X21S1+...
                                                   G_dot(4,2)*S2+...
                                                   G(3,2)*(-v2_relx*X21S1);
                          ];
 
                Sx = [obj.gearRatio*Sr1x    z;
                      z                     obj.gearRatio*Sr2x;
                      G(3,1)*S1x            G(3,2)*S1x;
                      G(3,1)*X21S1x+...
                      G(4,1)*S2x            G(3,2)*X21S1x+G(4,2)*S2x;
                     ];

                Sd = S_ring + Sx;

            end
        end

        function y_dot = get_minimal_coordinate_velocity(obj, ql, ql_dot)
            % TODO (@nicholasadr): Hack to get minimal coordinate in ID and
            %                      FDab.
            y_dot_fn = casadi.external(...
                            strcat(obj.type,'_IK_vel'),...
                            strcat('./models/casadi_fn/gen_',...
                                   obj.type,'_IK_dql_to_dy.so'));
            y_dot = full( y_dot_fn(ql,ql_dot) );
        end
    end
end