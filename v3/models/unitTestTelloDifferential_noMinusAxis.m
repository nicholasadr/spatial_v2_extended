clear
%% Test Tello Differential
disp('Tello Differential Unit Test No Minus Axis')

model = tello_model();

% Random configuration and velocity
% - ql: dependent maximal coordinate position (link joint position)
% - ql_dot: dependent maximal coordinate velocity (link joint velocity)
% - y_ddot: minimal coordinate acceleration (post-gearbox rotor acceleration)
% note: ql, ql_dot are also the minimal coordinates for cluster of index
%       2,5
fb_q = [1 0 0 0 0 0 0]'; % [a, bi, cj, dk, x, y, z] of quaternion rotation (a+bi+cj+dk) and x,y,z position
fb_qd = [0 0 0 0 0 0]';
fb_qdd = [0 0 0 0 0 0]';
ql   = [fb_q; rand(model.NQ-7,1)];
ql   = normalizeConfVec(model, ql);
ql_dot  = [fb_qd; rand(model.NV-6,1)];
y_ddot = [fb_qdd; rand(model.NV-6,1)];

% Uncomment below to check fixed configuration
% Spanning tree (position) config for differential joints:
% left hip differential: 0.9063 0.8797 8.7312 -0.6608
% left knee ankle differential: 0.8178 0.2607 3.4109 6.4027
% right hip differential: 0.0225 0.4253 3.7552 3.4770
% right knee ankle differential: 0.3127 0.1615 0.9567 2.7957
ql = [1. 0. 0. 0. 0. 0. 0. 0.1078 0.9063 0.8797 0.8178 0.2607 0.5944 0.0225 0.4253 0.3127 0.1615]';
ql_dot = [0. 0. 0. 0. 0. 0. 0.1788 0.4229 0.0942 0.5985 0.4709 0.6959 0.6999 0.6385 0.0336 0.0688]';
y_ddot = [0. 0. 0. 0. 0. 0. 0.3196 0.5309 0.6544 0.4076 0.8200 0.7184 0.9686 0.5313 0.3251 0.1056]';
% (for cpp) q_dot = [0. 0. 0. 0. 0. 0. 0.1788 -0.6479 -0.4486 0.1383 1.0587 0.6959 1.4263 -0.0333 -0.0325 0.0997]';
       
% Calculate dynamics quantities
[tau, out] = ID(model, ql, ql_dot ,y_ddot);  % Inverse dynamics
q_ddot_ABA = FDab(model, ql, ql_dot, tau);   % Forward dynamics