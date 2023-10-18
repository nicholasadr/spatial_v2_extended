clear
%% Test Tello Differential
disp('Tello Differential Unit Test No Minus Axis')

import casadi.*

model = RBD_model();
% uncomment below to switch off gravity
%disp('gravity: off')
%model.gravity = [0;0;0];

model.NB = 7;
model.parent = [0 1 2 3 1 5 6];

gearRatio = 6;

% import tello inertial pararmeters from solidworks-generated data files
inertial = fcn_import_params_solidworks('tello_params/inertial_*.txt');

model.joint{1} = floatingBaseJoint();
model.Xtree{1} = eye(6);
model.I{1} = mcI(inertial.mass{1},...
                 inertial.CoM{1},...
                 inertial.Ic{1});

rotor_mass = 0.07;
rotor_CoM = [0 0 0]';
rotor_Ixx = 2.5984e-05;
rotor_Iyy = rotor_Ixx;
rotor_Izz = 5.1512e-05;
I_rotor = [rotor_Ixx rotor_Iyy rotor_Izz]' .* eye(3);
rotorInertia = mcI(rotor_mass, rotor_CoM, I_rotor);

% joint orientations
R_top = eye(3);
R_down = [1, 0, 0;...
          0, -1, 0;...
          0, 0, -1];
R_left = [-1, 0, 0;...
           0, 0, 1;...
           0, 1, 0];
R_right = [1, 0, 0;...
           0, 0, 1;...
           0, -1, 0];

% 1st cluster: left hip clamp
model.joint{2} = revoluteJointWithRotorTello();
model.joint{2}.gearRatio = gearRatio;
model.joint{2}.rotorAxis = 'Rz';
model.joint{2}.jointAxis = 'Rz';
model.Xtree{2}(1:6,:) = plux(R_down, [0, 126, -26]*1E-3); % rotor
model.Xtree{2}(7:12,:) = plux(R_top, [0, 126, -87]*1E-3); % link
model.I{2}(1:6,1:6) = rotorInertia;
model.I{2}(7:12,7:12) = mcI(inertial.mass{2},...
                            inertial.CoM{2},...
                            inertial.Ic{2}); % hip clamp link

% 2nd cluster: left hip differential
model.joint{3} = telloDifferential();
model.joint{3}.type = 'thd';
model.joint{3}.rotorAxis{1} = 'Rz';
model.joint{3}.rotorAxis{2} = 'Rz';
model.joint{3}.jointAxis{1} = 'Rx';
model.joint{3}.jointAxis{2} = 'Ry';
model.joint{3}.XtreeInternal = plux(R_top, [0,0,0]*1E-3); % from first link (gimbal) to second link (thigh)
model.Xtree{3}(1:6,:) = plux(R_left, [0, 40, 0]*1E-3); % from predecessor to first rotor
model.Xtree{3}(7:12,:) = plux(R_right, [0, -40, 0]*1E-3); % from predecessor to second rotor
model.Xtree{3}(13:18,:) = plux(R_top, [0, 0, -142.5]*1E-3); % from predecessor to first link (gimbal)
model.I{3}(1:6,1:6) = rotorInertia; % 1st rotor
model.I{3}(7:12,7:12) = rotorInertia; % 2nd rotor 
model.I{3}(13:18,13:18) = mcI(inertial.mass{3},...
                              inertial.CoM{3},...
                              inertial.Ic{3}); % 1st link (gimbal)
model.I{3}(19:24,19:24) = mcI(inertial.mass{4},...
                              inertial.CoM{4},...
                              inertial.Ic{4}); % 2nd link (thigh)

% 3rd cluster: left knee-ankle differential
model.joint{4} = telloDifferential();
model.joint{4}.type = 'tkad';
model.joint{4}.rotorAxis{1} = 'Rz';
model.joint{4}.rotorAxis{2} = 'Rz';
model.joint{4}.jointAxis{1} = 'Ry';
model.joint{4}.jointAxis{2} = 'Ry';
model.joint{4}.XtreeInternal = plux(R_top, [0, 0, -260]*1E-3); % from first link (shank) to second link (foot)
model.Xtree{4}(1:6,:) = plux(R_right, [0, 26.55, 0]*1E-3); % from predecessor to first rotor
model.Xtree{4}(7:12,:) = plux(R_left, [0, -26.55, 0]*1E-3); % from predecessor to second rotor
model.Xtree{4}(13:18,:) = plux(R_top, [0, 0, -226.8]*1E-3); % from predecessor to first link (shank)
model.I{4}(1:6,1:6) = rotorInertia; % 1st rotor 
model.I{4}(7:12,7:12) = rotorInertia; % 2nd rotor 
model.I{4}(13:18,13:18) = mcI(inertial.mass{5},...
                              inertial.CoM{5},...
                              inertial.Ic{5}); % 1st link (shank)
model.I{4}(19:24,19:24) = mcI(inertial.mass{6},...
                              inertial.CoM{6},...
                              inertial.Ic{6}); % 2nd link (foot)

% 4th cluster: right hip clamp
model.joint{5} = revoluteJointWithRotorTello();
model.joint{5}.gearRatio = gearRatio;
model.joint{5}.rotorAxis = 'Rz';
model.joint{5}.jointAxis = 'Rz';
model.Xtree{5}(1:6,:) = plux(R_down, [0, -126, -26]*1E-3); % rotor
model.Xtree{5}(7:12,:) = plux(R_top, [0, -126, -87]*1E-3); % link
model.I{5}(1:6,1:6) = rotorInertia; 
model.I{5}(7:12,7:12) = mcI(inertial.mass{2},...
                            inertial.CoM{2},...
                            inertial.Ic{2}); % hip clamp link

% 5th cluster: right hip differential
model.joint{6} = telloDifferential();
model.joint{6}.type = 'thd';
model.joint{6}.rotorAxis{1} = 'Rz';
model.joint{6}.rotorAxis{2} = 'Rz';
model.joint{6}.jointAxis{1} = 'Rx';
model.joint{6}.jointAxis{2} = 'Ry';
model.joint{6}.XtreeInternal = plux(R_top, [0,0,0]*1E-3); % from first link (gimbal) to second link (thigh)
model.Xtree{6}(1:6,:) = plux(R_left, [0, 40, 0]*1E-3); % from predecessor to first rotor
model.Xtree{6}(7:12,:) = plux(R_right, [0, -40, 0]*1E-3); % from predecessor to second rotor
model.Xtree{6}(13:18,:) = plux(R_top, [0, 0, -142.5]*1E-3); % from predecessor to first link (gimbal)
model.I{6}(1:6,1:6) = rotorInertia; % 1st rotor 
model.I{6}(7:12,7:12) = rotorInertia; % 2nd rotor 
model.I{6}(13:18,13:18) = mcI(inertial.mass{3},...
                              inertial.CoM{3},...
                              inertial.Ic{3}); % 1st link (gimbal)
model.I{6}(19:24,19:24) = mcI(inertial.mass{4},...
                              inertial.CoM{4},...
                              inertial.Ic{4}); % 2nd link (thigh)

% 6th cluster: right knee-ankle differential
model.joint{7} = telloDifferential();
model.joint{7}.type = 'tkad';
model.joint{7}.rotorAxis{1} = 'Rz';
model.joint{7}.rotorAxis{2} = 'Rz';
model.joint{7}.jointAxis{1} = 'Ry';
model.joint{7}.jointAxis{2} = 'Ry';
model.joint{7}.XtreeInternal = plux(R_top, [0, 0, -260]*1E-3); % from first link (shank) to second link (foot)
model.Xtree{7}(1:6,:) = plux(R_right, [0, 26.55, 0]*1E-3); % from predecessor to first rotor
model.Xtree{7}(7:12,:) = plux(R_left, [0, -26.55, 0]*1E-3); % from predecessor to second rotor
model.Xtree{7}(13:18,:) = plux(R_top, [0, 0, -226.8]*1E-3); % from predecessor to first link (shank)
model.I{7}(1:6,1:6) = rotorInertia; % 1st rotor 
model.I{7}(7:12,7:12) = rotorInertia; % 2nd rotor 
model.I{7}(13:18,13:18) = mcI(inertial.mass{5},...
                              inertial.CoM{5},...
                              inertial.Ic{5}); % 1st link (shank)
model.I{7}(19:24,19:24) = mcI(inertial.mass{6},...
                              inertial.CoM{6},...
                              inertial.Ic{6}); % 2nd link (foot)

model = model.postProcessModel();

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