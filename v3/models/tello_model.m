function model = tello_model(gravity)

import casadi.*

% NOTE: TODO

if nargin == 1
    model.gravity = [0; 0; gravity];
end

model = RBD_model();

model.NB = 7;
model.parent = [0 1 2 3 1 5 6];
mode.appearance.base = {'tiles', [-10 10; -10 10; 0 0], 0.1,'colour',[0.047 0.137 0.25],'box',[-0.07 0 -0.01; 0.15 0.05 0.01]};

gearRatio = 6;

% import tello inertial parameters from solidworks-generated data files
inertial = fcn_import_params_solidworks('tello_params/inertial_*.txt');

% 1st cluster: torso
model.joint{1} = floatingBaseJoint();
model.Xtree{1} = eye(6);
model.I{1} = mcI(inertial.mass{1},...
                 inertial.CoM{1},...
                 inertial.Ic{1});
model.appearance.body{1} = {{'3dobj',"obj/tello/torso.obj"}};

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

% 2nd cluster: left hip clamp
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
model.appearance.body{2} = {{},{'3dobj',"obj/tello/hipclamp.obj"}};

% 3rd cluster: left hip differential
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
model.appearance.body{3} = {{},{},{'3dobj',"obj/tello/gimbal.obj"},...
                            {'3dobj',"obj/tello/thigh.obj"}};

% 4th cluster: left knee-ankle differential
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
model.appearance.body{4} = {{},{},{'3dobj',"obj/tello/shin.obj"},...
                            {'3dobj',"obj/tello/foot.obj"}};
% 5th cluster: right hip clamp
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
model.appearance.body{5} = {{},{'3dobj',"obj/tello/hipclamp.obj"}};

% 6th cluster: right hip differential
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
model.appearance.body{6} = {{},{},{'3dobj',"obj/tello/gimbal.obj"},...
                           {'3dobj',"obj/tello/thigh.obj"}};

% 7th cluster: right knee-ankle differential
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
model.appearance.body{7} = {{},{},{'3dobj',"obj/tello/shin.obj"},...
                            {'3dobj',"obj/tello/foot.obj"}};

model.camera.direction = [1 0 0];
model.camera.up = [0 0 1];

model = model.postProcessModel();

end
