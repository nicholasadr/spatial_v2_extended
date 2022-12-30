function tello_visualize_2()
% todo - optional: rotor visualization 

%  visualize structural bodies only
obj_filename{1} = '.\vis\torso.obj';

obj_filename{2} = '.\vis\hipclamp.obj';
obj_filename{3} = '.\vis\gimbal.obj';
obj_filename{4} = '.\vis\thigh.obj';
obj_filename{5} = '.\vis\shin.obj';
obj_filename{6} = '.\vis\foot.obj';

obj_filename{12} = '.\vis\hipclamp.obj';
obj_filename{13} = '.\vis\gimbal.obj';
obj_filename{14} = '.\vis\thigh.obj';
obj_filename{15} = '.\vis\shin.obj';
obj_filename{16} = '.\vis\foot.obj';

obj_filename{22} = '.\vis\spatial_4bar_top.obj';
obj_filename{23} = '.\vis\spatial_4bar_top.obj';
obj_filename{24} = '.\vis\carrier.obj';
obj_filename{25} = '.\vis\hamstring.obj';

obj_filename{26} = '.\vis\spatial_4bar_top.obj';
obj_filename{27} = '.\vis\spatial_4bar_top.obj';
obj_filename{28} = '.\vis\carrier.obj';
obj_filename{29} = '.\vis\hamstring.obj';

% read mesh files as vertices and faces
for ii = [1, 2:6, 12:16, 22:29]
    obj{ii} = readObj(obj_filename{ii});
end

% draw robot in this configuration
% q = zeros(29, 1);
% q = get_config([0;0;-2*pi/6;3.5*pi/6;0]);
% q = get_config([0;0;0;pi/6;0]);
% q = get_config([0;0;0;0;0])
q = get_config([0;0;0;pi/3;0]);

[model]     = tello_w_rotors_w_parallel();
NB = model.NB;
p  = model.parent;

color_ = parula(30);
color = color_(1:NB,:);
%     color_rot = color_(NB+1:end,:);

hold on;
body_index = [1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 28, 29];
for i = body_index
    [ XJ, ~]     = jcalc( model.jtype{i}, q(i) );
    %       [ XJrot, ~]  = jcalc( model.jtype{i}, q(i)*model.gr{i} );
    %       [ XJ_CAD, ~] = jcalc( model.jtype{i}, q_CAD(i) );

    %       if i == 2
    %          % Account for offset CAD origin for thigh link
    %          XJ_CAD = XJ_CAD * [eye(3,3),zeros(3,3);skew([0;-0.045;0]),eye(3,3)];
    %       end

    Xup{i} = XJ * model.Xtree{i};
    %       Xrot{i} = XJrot * model.Xrotor{i};
    %       Xup_CAD{i} = XJ_CAD * model.Xtree{i};

    if model.parent(i) == 0
        X0{i}     = Xup{i};
        %         X0rot{i}  = Xrot{i};
        %         X0_CAD{i} = Xup_CAD{i};

    else
        X0{i}     = Xup{i} *  X0{p(i)};
        %         X0rot{i}  = Xrot{i} * X0{p(i)};
        %         X0_CAD{i} = Xup_CAD{i} *  X0{p(i)};

    end

    T    = inv( AdjointToSE3( X0{i} ) );
    %       Trot = inv( AdjointToSE3( X0rot{i} ) );
    %       T = inv( AdjointToSE3( X0_CAD{i} ) );

    if i>0
%         if (i == 24) || (i == 28)
%             R = ry(-40/180*pi);
%             obj{i}.v = (R*obj{i}.v')';
%         end
        V = (T(1:3,1:3)*obj{i}.v' + T(1:3,4) * ones(1,size(obj{i}.v,1)))';
        F = obj{i}.f.v;
        patch(  'vertices', V, ...
            'faces', F, ...
            'LineStyle','-', ...
            'EdgeAlpha', .4, ...
            'EdgeColor', color(i, :)*.95, ...
            'FaceAlpha',.2, ...
            'FaceColor', color(i, :));
    end
end
axis equal;
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, .5])
grid on
view([-140 26]);
draw_SE3(eye(4,4));
end


function qm = get_config(q)
qm = zeros(29, 1);
qm(1) = 0;
qm(2:6) = q;
ik = fcn_ik_q_2_p(q);
qm(7:11) = ik;
qm(12:16) = q;
qm(17:21) = fcn_ik_q_2_p(q);
qm(22) = ik(2);
qm(23) = ik(3);
qm(26) = ik(2);
qm(27) = ik(3);

qm(24) = q(3) + q(4);
qm(25) = -q(4);
qm(28) = q(3) + q(4);
qm(29) = -q(4);
end

