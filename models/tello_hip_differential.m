function robot = tello_hip_differential(qd, y_dot, y_ddot)
% qd: 2x1 vector of dependent joint coordinates
% (hip roll and hip pitch joints)
% y_dot: 2x1 vector of independent joint velocities (actuators)
% y_ddot: 2x1 vector of independent joint accelerations (actuators)
% Example to run:
% >> qd = [-0.5059 -0.1865]';
% >> y_dot = [-0.3632 1.2884]';
% >> y_ddot = [0 0]';
% >> hip = tello_hip_differential(qd, y_dot, y_ddot);
import casadi.*

% NOTE 1: Inertial parameters are stored in 
% 'inertial_params_X_partname.txt'. These text files are manually generated 
% from solidworks. 
% NOTE 2: Following cases are added to jcalc
%     case '-Rx'
%         Xj = rotx(q);
%         S = [-1;0;0;0;0;0];
%     case '-Ry'				% revolute Y axis
%         Xj = roty(-q);
%         S = [0;-1;0;0;0;0];
%     case '-Rz'			% revolute Z axis
%         Xj = rotz(-q);
%         S = [0;0;-1;0;0;0];
% NOTE 3: Body index:
%   0: hip clamp
%   1: rotor 2
%   2: rotor 3
%   3: gimbal
%   4: thigh

% number of bodies
NB = 4; 

% tree structure
parent = [0, 0, 0, 3];

% joint types
jtype = { 'Rz', 'Rz', 'Rx', 'Ry' };

% joint orientations
R2 = [  -1, 0, 0; ...
        0, 0, 1; ...
        0, 1, 0];

% joint tree transformations
Xtree = { 
         plux(R2, [0, 40, 0]*1E-3), ...
         plux(R2, [0, -40, 0]*1E-3), ...
         plux(eye(3), [0, 0, -142.5]*1E-3), ...
         plux(eye(3), [0,0,0]*1E-3), ...
         };

% spatial inertia
I{1} = mcI( 1, ...
            [0 0 0]', ...
            [0 0 0; 0 0 0; 0 0 0.003]);
I{2} = mcI( 1, ...
            [0 0 0]', ...
            [0 0 0; 0 0 0; 0 0 0.003]);
inertial_3 = fcn_import_params_solidworks('./models/tello_params/inertial_params_3_gimbal.txt');
inertial_4 = fcn_import_params_solidworks('./models/tello_params/inertial_params_4_thigh.txt');
I{3} = mcI( inertial_3.mass{1}, ...
            inertial_3.CoM{1}', ...
            inertial_3.Ic{1});
I{4} = mcI( inertial_4.mass{1}, ...
            inertial_4.CoM{1}', ...
            inertial_4.Ic{1});

% inverse kinematics of hip differential: qd -> y
y_fn = external('IK_dependent_state_to_y', './models/casadi_fn/y_gen.so');
y = full(y_fn(qd));
q = [y(1) y(2) qd(1) qd(2)]';

for i = 1:NB
    [ XJ{i}, S_single_joint{i} ] = jcalc( jtype{i}, q(i) );
    Xup{i} = XJ{i} * Xtree{i};
end

% jacobian that maps independent velocities (actuators) to dependent
% velocities (joints): y_dot -> qd_dot
J_q_hip = (zeros(2,2));
J_q_hip(1,1) = ...
    (49*sin(qd(1)))/5000 - (399*cos(qd(1)))/20000 - (7*cos(qd(1))*sin(y(1)))/625 + ...
    (7*cos(qd(1))*sin(qd(2)))/625 + (57*sin(qd(1))*sin(qd(2)))/2500 + ...
    (8*sin(y(1))*sin(qd(1))*sin(qd(2)))/625;
J_q_hip(1,2) = ...
    (8*cos(y(1))*sin(qd(2)))/625 - ...
    (57*cos(qd(1))*cos(qd(2)))/2500 + (7*cos(qd(2))*sin(qd(1)))/625 - ...
    (8*cos(qd(1))*cos(qd(2))*sin(y(1)))/625;
J_q_hip(2,1) = ...
    (399*cos(qd(1)))/20000 + (49*sin(qd(1)))/5000 + (7*cos(qd(1))*sin(y(2)))/625 - ...
    (7*cos(qd(1))*sin(qd(2)))/625 + (57*sin(qd(1))*sin(qd(2)))/2500 + ...
    (8*sin(y(2))*sin(qd(1))*sin(qd(2)))/625;
J_q_hip(2,2) = ...
    (8*cos(y(2))*sin(qd(2)))/625 - (57*cos(qd(1))*cos(qd(2)))/2500 - ...
    (7*cos(qd(2))*sin(qd(1)))/625 - (8*cos(qd(1))*cos(qd(2))*sin(y(2)))/625;
 
J_p_hip = (zeros(2,2));
J_p_hip(1,1) = ...
    (57*cos(y(1)))/2500 - (7*cos(y(1))*sin(qd(1)))/625 + (8*cos(qd(2))*sin(y(1)))/625 - ...
    (8*cos(y(1))*cos(qd(1))*sin(qd(2)))/625; ...
J_p_hip(2,2) = ...
    (57*cos(y(2)))/2500 + (7*cos(y(2))*sin(qd(1)))/625 + ...
    (8*cos(qd(2))*sin(y(2)))/625 - (8*cos(y(2))*cos(qd(1))*sin(qd(2)))/625;
 
J_dp_2_dq_hip = -J_q_hip\J_p_hip;

G = [1 0; 0 1; J_dp_2_dq_hip];

a_grav = [0 0 0 0 0 -9.81]';
q_dot = G * y_dot;
g_fn = external('g_gen', './models/casadi_fn/g_gen.so');
g = full(g_fn(y,qd,[q_dot(1);q_dot(2)],[q_dot(3);q_dot(4)]));
q_ddot = G * y_ddot + g;

% RNEA (Rigid Body Tree Model)
for i = 1:NB
  %[ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ{i} = S_single_joint{i}*q_dot(i);
  Xup{i} = XJ{i} * Xtree{i};
  if parent(i) == 0
    v{i} = vJ{i};
    a{i} = Xup{i}*(-a_grav) + S_single_joint{i}*q_ddot(i); %S_ring is zero, v*vJ is zero
  else
    v{i} = Xup{i}*v{parent(i)} + vJ{i};
    a{i} = Xup{i}*a{parent(i)} + S_single_joint{i}*q_ddot(i) + crm(v{i})*vJ{i}; %S_ring is zero
  end
  f{i} = I{i}*a{i} + crf(v{i})*I{i}*v{i};
end

for i = NB:-1:1
  tau(i,1) = S_single_joint{i}' * f{i};
  if parent(i) ~= 0
    f{parent(i)} = f{parent(i)} + Xup{i}'*f{i};
  end
end

% spatial inertia of cluster
I_ = zeros(NB*6,NB*6);
I_(1:6,1:6) = I{1};
I_(7:12,7:12) = I{2};
I_(13:18,13:18) = I{3};
I_(19:24,19:24) = I{4};

parent_ = 0;

% Cluster Xup
Xup_ = zeros(NB*6,6);
Xup_(1:6,:) = Xup{1};
Xup_(7:12,:) = Xup{2};
Xup_(13:18,:) = Xup{3};
Xup_(19:24,:) = Xup{4}*Xup{3};

% Calculate S
% - X_span
X_span = zeros(NB*6,NB*6);
X_span(1:6,1:6) = eye(6);
X_span(7:12,7:12) = eye(6);
X_span(13:18,13:18) = eye(6);
X_span(19:24,13:18) = Xup{4};
X_span(19:24,19:24) = eye(6);
% - S_span
S_span = zeros(NB*6,4);
S_span(1:6,1) = S_single_joint{1};
S_span(7:12,2) = S_single_joint{2};
S_span(13:18,3) = S_single_joint{3};
S_span(19:24,4) = S_single_joint{4};
% - S
S = X_span * S_span * G;

% calculate S_ring
% - calculate a,b,c,d
abcd = J_dp_2_dq_hip;
abcd_ = abcd';
abcd_ = abcd_(:)';
tmp = num2cell(abcd_);
[A B C D] = deal(tmp{:});
% - calculate a_dot, b_dot, c_dot, d_dot
kikd_fn = external('kikd_gen', './models/casadi_fn/kikd_gen.so');
[Ki,Kd,Ki_dot,Kd_dot] = kikd_fn([y(1);y(2)],[qd(1);qd(2)],[q_dot(1);q_dot(2)],[q_dot(3);q_dot(4)]);
abcd_dot = -Kd\Ki_dot + (Kd\Kd_dot)*(Kd\Ki);
abcd_dot = full(abcd_dot); % convert from casadi.DM to array
abcd_dot_ = abcd_dot';
abcd_dot_ = abcd_dot_(:)';
tmp = num2cell(abcd_dot_);
[A_dot B_dot C_dot D_dot] = deal(tmp{:});
% S_ring
S_ring = zeros(NB*6,2);
S_ring(13:18,1) = A_dot * S_single_joint{3};
S_ring(13:18,2) = B_dot * S_single_joint{3};
S_ring(19:24,1) = A_dot * Xup{4} * S_single_joint{3} +...
                  C_dot * S_single_joint{4} +...
                  A * -crm(vJ{4}) * Xup{4} * S_single_joint{3};
S_ring(19:24,2) = B_dot * Xup{4} * S_single_joint{3} +...
                  D_dot * S_single_joint{4} +...
                  B * -crm(vJ{4}) * Xup{4} * S_single_joint{3};

% RNEA (Cluster Tree Model)
%[ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
vJ_ = S*y_dot;
%Xup{i} = XJ{i} * Xtree{i};
if parent_ == 0
  v_ = vJ_;
  a_ = Xup_*(-a_grav) + S*y_ddot + S_ring*y_dot; %v*vJ is zero
%else
%  v_ = Xup_*v{parent(i)} + vJ{i};
%  a_ = Xup_*a{parent(i)} + S_single_joint{i}*q_ddot(i) + crm(v{i})*vJ{i};
end
crfv = zeros(NB*6,NB*6);
crfv(1:6,1:6) = crf(v_(1:6));
crfv(7:12,7:12) = crf(v_(7:12));
crfv(13:18,13:18) = crf(v_(13:18));
crfv(19:24,19:24) = crf(v_(19:24));
%f_ = I_*a_ + crf(v_)*I_*v_;
f_ = I_*a_ + crfv*I_*v_;
tau_ = S' * f_;

%% debug
%G3 = [A B];
%G3_dot = [A_dot B_dot];
%a1_ = S_single_joint{3}*G3_dot*y_dot; % [-0.2114; 0;0;0;0;0]
%a2_ = S_ring(13:18,:)*y_dot; % [-0.2114; 0;0;0;0;0]
%a3_ = G3_dot * y_dot; % -0.2114
%a4_ = g(3); % 0.1075

robot.NB = NB;
robot.parent = parent;
robot.jtype = jtype;
robot.Xtree = Xtree;
robot.XJ = XJ;
robot.Xup = Xup;
robot.S_single_joint = S_single_joint;
robot.G = G;
robot.S = S;
robot.S_ring = S_ring;
robot.vJ = vJ;
robot.v = v;
robot.a = a;
robot.f = f;
robot.tau = tau;
robot.vJ_ = vJ_;
robot.v_ = v_;
robot.a_ = a_;
robot.f_ = f_;
robot.tau_ = tau_;
robot.q = q;
robot.q_dot = q_dot;
robot.q_ddot = q_ddot;
end