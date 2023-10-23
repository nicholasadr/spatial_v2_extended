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
    [q] = confVecToCell(model,q);
end

for i = 1:model.NB
  [Xup{i}, S{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i});
  if model.parent(i) == 0
      X0{i} = Xup{i};
  else
      X0{i} = Xup{i} * X0{model.parent(i)};
  end
end

for i = 1:size(model.bodies_with_obj_file_inds,2)
    b_jk = model.bodies_with_obj_file_inds(:,i);
    cluster_idx = b_jk(1); % index of cluster
    in_cluster_idx = b_jk(2); % index of body in cluster
    if i == 1
        X0_body = X0{i};
    else
        row_inds = 6*(in_cluster_idx-1)+1 : 6*in_cluster_idx;
        X0_body = X0{cluster_idx}(row_inds, :);
    end

    T = inv(pluho(X0_body));
    obj = readObj(model.obj{cluster_idx}(in_cluster_idx));
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

axis equal;
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, .5])
view([-140 26]);
grid on
draw_SE3(eye(4,4));
end

