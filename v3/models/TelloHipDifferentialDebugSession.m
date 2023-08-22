qd = [-0.5059 -0.1865]';
qd_dot = [-0.9651 0.4352]';
qd_ddot = [   -0.2114 0.1344]';

hipd_dep = tello_hip_differential_from_dependent(qd, qd_dot, qd_ddot);

a1_cluster_dep = hipd_dep.a_(13:18,:)

a1_tree_dep = hipd_dep.a{3}


G_dep = hipd_dep.G
g_dep = hipd_dep.g

Jinv = G_dep(3:4,:);
J = inv(Jinv);


y_dot = [-0.3632 1.2884]';
y_ddot = [0 0]';
hipd_indep = tello_hip_differential(qd, y_dot, y_ddot);

a1_cluster_indep = hipd_indep.a_(13:18)
a1_tree_indep    = hipd_indep.a{3}

G_indep = hipd_indep.G
g_indep = hipd_indep.g


qd_ddot_1 = G_indep*y_ddot + g_indep
qd_ddot_2 = [y_ddot ; qd_ddot]

% J*qd_dot = Y_dot
%
J = full(f_jac(qd))
Jdot = full(f_jac_dot(qd,qd_dot))
y_ddot2 = J*qd_ddot + Jdot*qd_dot
y_ddot1 = full( f3(qd,qd_dot, qd_ddot) )

% 
