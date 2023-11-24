function tello_visualize_grbda()

model = tello_model();

no_of_unique_color = 6;
color_ = parula(no_of_unique_color);
color = color_(1:no_of_unique_color,:);
color = [color; color(2:end,:)];

% TODO: create a function for the forward kinematics below

a_grav = model.getGravity();

q_fb = [1 0 0 0 0 0 0]';
q = [q_fb; zeros(model.NQ,1)];

if ~iscell(q)
    [q] = model.confVecToCell(q);
end

for i = 1:model.NB
  [Xup{i}, S{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i});
  if model.parent(i) == 0
      X0{i} = Xup{i};
  else
      X0{i} = Xup{i} * X0{model.parent(i)};
  end
end

for cluster_idx = 1:model.NB
    if i == 1
      X0_body = X0{i};
    else
      in_cluster_idx_list = find(cell2mat(model.obj_bool_mask{cluster_idx})); % indices of bodies in cluster with 3d obj
      for in_cluster_idx = in_cluster_idx_list'
        row_inds = 6*(in_cluster_idx-1)+1 : 6*in_cluster_idx;
        X0_body = X0{cluster_idx}(row_inds, :);

        T = inv(pluho(X0_body));
        obj = readObj(model.appearance.body{cluster_idx}{in_cluster_idx}{end});
        V = (T(1:3,1:3)*obj.v' + T(1:3,4) * ones(1,size(obj.v,1)))';
        F = obj.f.v;
        patch('vertices', V, ...
          'faces', F, ...
          'LineStyle','-', ...
          'EdgeAlpha', .4, ...
          'EdgeColor', color(i, :)*.95, ...
          'FaceAlpha',.2, ...
          'FaceColor', color(i, :));
      end
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

