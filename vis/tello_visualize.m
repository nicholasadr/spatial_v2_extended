function tello_visualize()
% todo: rotor visualization 

% only visualize structural bodies
obj_filename{1} = 'torso.obj';

obj_filename{2} = 'hipclamp.obj';
obj_filename{3} = 'gimbal.obj';
obj_filename{4} = 'thigh.obj';
obj_filename{5} = 'shin.obj';
obj_filename{6} = 'foot.obj';

obj_filename{12} = 'hipclamp.obj';
obj_filename{13} = 'gimbal.obj';
obj_filename{14} = 'thigh.obj';
obj_filename{15} = 'shin.obj';
obj_filename{16} = 'foot.obj';

% read mesh files as vertices and faces
for ii = [1, 2:6, 12:16]
    obj{ii} = readObj(obj_filename{ii});
end

% draw robot in this configuration
q = zeros(21, 1);

[model]     = tello_w_rotors_wo_parallel();
NB = model.NB;
p  = model.parent;

color_ = parula(25);
color = color_(1:NB,:);
%     color_rot = color_(NB+1:end,:);

hold on;
body_index = [1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 16];
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
view([-140 26]);
grid on
draw_SE3(eye(4,4));
end