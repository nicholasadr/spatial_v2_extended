function model = mit_humanoid_model(gravity)
    % x forward, y left, z up

import casadi.*

% NOTE: TODO

if nargin == 1
    model.gravity = [0; 0; gravity];
end

model = RBD_model();

model.NB = 17;
%torso
%leg: [hiprz, hiprz_rotor], [hiprx, hiprx_rotor], [hipry, hipry_rotor], [knee, knee_rotor, ankle, ankle_rotor]
% arm: [shoulderRy, shoulderRy_rotor], [shoulderRx, shoulderRx_rotor], [shoulderRz, shoulderRz_rotor], [elbow, elbow_rotor]
model.parent = [0 1 2 3 4 1 6 7 8 1 10 11 12 1 14 15 16];
%model.appearance.base = {'tiles', [-10 10; -10 10; 0 0], 0.1,'colour',[0.047 0.137 0.25],'box',[-0.07 0 -0.01; 0.15 0.05 0.01]};

% Link rotational inertia
torsoRotInertia = [0.172699, 0.001419, 0.004023;...
                   0.001419, 0.105949, -0.001672;...
                   0.004023, -0.001672, 0.091906];
hipRzRotInertia = [0.0015373, 0.0000011, 0.0005578;...
                   0.0000011, 0.0014252, 0.0000024;...
                   0.0005578, 0.0000024, 0.0012028];
hipRxRotInertia = [0.0017535, -0.0000063, -0.000080;...
                   -0.0000063, 0.003338, -0.000013;...
                   -0.000080, -0.000013, 0.0019927];
hipRyRotInertia = [0.0243761, 0.0000996, 0.0006548;...
                   0.0000996, 0.0259015, 0.0026713;...
                   0.0006548, 0.0026713, 0.0038929];
kneeRotInertia = [0.003051, 0.000000, 0.0000873;...
                  0.000000, 0.003033, 0.0000393;...
                  0.0000873, 0.0000393, 0.0002529];
ankleRotInertia = [0.0000842, 0.000000, -0.0000488;...
                   0.000000, 0.0007959, -0.000000;...
                   -0.0000488, -0.000000, 0.0007681];
shoulderRyRotInertia = [0.0013678, 0.0000266, 0.0000021;...
                	0.0000266, 0.0007392, -0.0000012;...
                	0.0000021, -0.0000012, 0.000884];
shoulderRxRotInertia = [0.0011524, 0.0000007, 0.0000396;...
                	0.0000007, 0.0011921, 0.0000014;...
                	0.0000396, 0.0000014, 0.0012386];
shoulderRzRotInertia = [0.0012713, 0.000001, -0.000008;...
                	0.000001, 0.0017477, -0.0000225;...
                	-0.000008, -0.0000225, 0.0008191];
elbowRotInertia = [0.001570, 0.0000002, 0.0000335;...
                   0.0000002, 0.0016167, 0.000003;...
                   0.0000335, 0.000003, 0.0000619];

% Rotor rotational inertia
largeRotorRotInertiaZ = [3.443e-4, 0, 0;...
                		  0, 3.443e-4, 0;...
                		  0, 0, 5.548e-4];
smallRotorRotInertiaZ = [1.084e-4, 0, 0;...
                		 0, 1.084e-4, 0;...
                		 0, 0, 1.6841e-4];
RY = roty(pi/2);
RY = RY(1:3,1:3);
RX = rotx(-pi/2);
RX = RX(1:3,1:3);

smallRotorRotInertiaX = RY' * smallRotorRotInertiaZ * RY;
smallRotorRotInertiaY = RX' * smallRotorRotInertiaZ * RX;
largeRotorRotInertiaY = RX' * largeRotorRotInertiaZ * RX;

% mass
torsoMass = 8.52;
hipRzMass = 0.84563;
hipRxMass = 1.208683;
hipRyMass = 2.64093;
kneeMass = 0.3543355;
ankleMass = 0.280951;
shoulderRyMass = 0.788506;
shoulderRxMass = 0.80125;
shoulderRzMass = 0.905588;
elbowMass = 0.34839;
rotorMass = 0.07; %TODO(@nicholasadr): added arbitrarily

% pitch for Xtree rotation
hipRzPitch = -0.174533;
hipRxPitch = 0.436332;
hipRyPitch = -(hipRxPitch + hipRzPitch);

% Xtree position
hipRzLocation = [-0.00565, -0.082, -0.05735]';
hipRxLocation = [-0.06435, 0.0, -.07499]';
hipRyLocation = [0.071, 0.0018375, 0.0]';
kneeLocation = [0.0, 0.0, -0.267]';
ankleLocation = [0.0, 0.0, -0.2785]';
shoulderRyLocation = [0.01346, -0.17608, 0.24657]';
shoulderRxLocation = [0.0, -0.0575, 0.0]';
shoulderRzLocation = [0.0, 0.0, -0.10250]';
elbowLocation = [0.0, 0.0, -0.1455]';

hipRzRotorLocation = [-0.00842837, -0.082, -0.041593]';
hipRxRotorLocation = [-0.0827, 0.0, -0.066436]';
hipRyRotorLocation = [0.071, 0.024, 0.0]';
kneeRotorLocation = [0.013, -0.0497, -0.0178]';
ankleRotorLocation = [0., 0., 0.]';
shoulderRyRotorLocation = [0.01346, -0.16, 0.24657]';
shoulderRxRotorLocation = [0, -0.0575, 0]';
shoulderRzRotorLocation = [0., 0., -0.1025]';
elbowRotorLocation = [0., 0.0325, -0.06]';

% CoM
torsoCOM = [0.009896, 0.004771, 0.100522]';

hipRzCOM = [-0.064842, -0.000036, -0.063090]';
hipRxCOM = [0.067232, -0.013018, 0.0001831]';
hipRyCOM = [0.0132054, 0.0269864, -0.096021]';
kneeCOM = [0.00528, 0.0014762, -0.13201]';
ankleCOM = [0.022623, 0.0, -0.012826]';

shoulderRyCOM = [0.009265, 0.052623, -0.0001249]';
shoulderRxCOM = [0.0006041, 0.0001221, -0.082361]';
shoulderRzCOM = [0.0001703, -0.016797, -0.060]';
elbowCOM = [-0.0059578, 0.000111, -0.0426735]';

smallRotorCOM = [0., 0., 0.]';
largeRotorCOM = [0., 0., 0.]';

% Gear ratio
hipRzGearRatio = 6.0;
hipRxGearRatio = 6.0;
hipRyGearRatio = 6.0;
kneeGearRatio = 6.0;
%kneeBeltRatios{2.0}; % {Thigh belt}
ankleGearRatio = 6.0;
%ankleBeltRatios{2.0, 1.0}; % {Thigh belt, Shank belt}
kneeBeltRatio = 2.0;
ankleBeltRatio = 1.0;

shoulderRxGearRatio = 6.0;
shoulderRzGearRatio = 6.0;
shoulderRyGearRatio = 6.0;
elbowGearRatio = 9.0;

% parameters for showmotion visualization
torso_length = 0.1;
torso_width = 0.164;
torso_height = 0.30392;
leg_rad = 0.02;
footToeLength = 0.1;
footHeelLength = 0.05;
footHeight = 0.041;
lowerArmLength = 0.27;

% 1st cluster: torso
model.joint{1} = floatingBaseJoint();
model.Xtree{1} = eye(6);
model.I{1} = mcI(torsoMass, torsoCOM, torsoRotInertia);
%model.appearance.body{1} = {{'3dobj','obj/mit_humanoid/torso.obj'}};
model.appearance.body{1} = {'colour',[232, 127, 35]./255,...
    'box', [-0.5*torso_length -0.75*torso_width -0.15*torso_height;...
    0.5*torso_length 0.75*torso_width 1.0*torso_height]};

idx = 1;
for legID = ["right","left"]

  % hipRz
  idx = idx + 1;
  model.joint{idx} = revoluteJointWithRotor();
  model.joint{idx}.gearRatio = hipRzGearRatio;
  model.joint{idx}.rotorAxis = 'z';
  model.joint{idx}.jointAxis = 'z';
  R_hipZ = roty(hipRzPitch);
  R_hipZ = R_hipZ(1:3,1:3);
  model.Xtree{idx}(1:6,:) = plux(R_hipZ, signedVector(hipRzLocation, legID));
  model.Xtree{idx}(7:12,:) = plux(R_hipZ, signedVector(hipRzRotorLocation, legID));
  model.I{idx}(1:6,1:6) = signedInertia(mcI(hipRzMass, hipRzCOM, hipRzRotInertia),...
					   legID);
  model.I{idx}(7:12,7:12) = signedInertia(mcI(rotorMass, smallRotorCOM, smallRotorRotInertiaZ),...
					     legID);
  %model.appearance.body{idx} = {{'3dobj','obj/mit_humanoid/hipRz.obj'},{}};
  cyl1 = [0 0 0.75*leg_rad;...
        0 0 -0.75*leg_rad];
  cyl2 = [0 0 0;...
          signedVector(hipRxLocation, legID)'];
  model.appearance.body{idx} = ...
        {{'colour',[18, 6, 186]./255,...
        'cyl', cyl1, 1.65*leg_rad,...
        'cyl', cyl2, 1.5*leg_rad}, {}};

  % hipRx
  idx = idx + 1;
  model.joint{idx} = revoluteJointWithRotor();
  model.joint{idx}.gearRatio = hipRxGearRatio;
  model.joint{idx}.rotorAxis = 'x';
  model.joint{idx}.jointAxis = 'x';
  R_hipX = roty(hipRxPitch);
  R_hipX = R_hipX(1:3,1:3);
  model.Xtree{idx}(1:6,:) = plux(R_hipX, signedVector(hipRxLocation, legID));
  model.Xtree{idx}(7:12,:) = plux(R_hipX, signedVector(hipRxRotorLocation, legID));
  model.I{idx}(1:6,1:6) = signedInertia(mcI(hipRxMass, hipRxCOM, hipRxRotInertia),...
					   legID);
  model.I{idx}(7:12,7:12) = signedInertia(mcI(rotorMass, smallRotorCOM, smallRotorRotInertiaX),...
					     legID);
  %model.appearance.body{idx} = {{'3dobj','obj/mit_humanoid/hipRx.obj'},{}};
  cyl1 = [0.5*leg_rad 0 0;...
        -1.5*leg_rad 0 0];
  cyl2 = [0 0 0;...
        signedVector(hipRyLocation, legID)'];
  model.appearance.body{idx} = ...
        {{'colour',[18, 6, 186]./255,...
        'cyl', cyl1, 2*leg_rad,...
        'cyl', cyl2, leg_rad}, {}};

  % hipRy
  idx = idx + 1;
  model.joint{idx} = revoluteJointWithRotor();
  model.joint{idx}.gearRatio = hipRyGearRatio;
  model.joint{idx}.rotorAxis = 'y';
  model.joint{idx}.jointAxis = 'y';
  R_hipY = roty(hipRyPitch);
  R_hipY = R_hipY(1:3,1:3);
  model.Xtree{idx}(1:6,:) = plux(R_hipY, signedVector(hipRyLocation, legID));
  model.Xtree{idx}(7:12,:) = plux(R_hipY, signedVector(hipRyRotorLocation, legID));
  model.I{idx}(1:6,1:6) = signedInertia(mcI(hipRyMass, hipRyCOM, hipRyRotInertia),...
					   legID);
  model.I{idx}(7:12,7:12) = signedInertia(mcI(rotorMass, largeRotorCOM, smallRotorRotInertiaY),...
					     legID);
  %model.appearance.body{idx} = {{'3dobj','obj/mit_humanoid/thigh.obj'},{}};
  cyl1 = [0 1.5*leg_rad 0;...
        0 -1.5*leg_rad 0];
  cyl2 = [0 0 0;...
        signedVector(kneeLocation, legID)'];
  model.appearance.body{idx} = ...
        {{'colour',[18, 6, 186]./255,...
        'cyl', cyl1, 3*leg_rad,...
        'cyl', cyl2, leg_rad}, {}};

  % knee-ankle
  idx = idx + 1;
  model.joint{idx} = revolutePairAbsoluteWithRotors();
  model.joint{idx}.gearRatio = {kneeGearRatio, ankleGearRatio};
  model.joint{idx}.beltRatio = {kneeBeltRatio; ankleBeltRatio}; %TODO(@nicholasadr): why is this not a scalar?
  model.joint{idx}.rotorAxis = {'y','y'};
  model.joint{idx}.jointAxis = {'y','y'};
  I3 = eye(3);
  model.joint{idx}.XtreeInternal = plux(I3, signedVector(ankleLocation, legID)); % from first link to second link;
  model.Xtree{idx}(1:6,:) = plux(I3, signedVector(kneeLocation, legID)); % from predecessor to first link
  model.Xtree{idx}(7:12,:) = plux(I3, signedVector(kneeRotorLocation, legID)); % from predecessor to first rotor
  model.Xtree{idx}(13:18,:) = plux(I3, signedVector(ankleRotorLocation, legID)); % from predecessor to second rotor
  model.I{idx}(1:6,1:6) = signedInertia(mcI(kneeMass, kneeCOM, kneeRotInertia), legID); % first link (knee)
  model.I{idx}(7:12,7:12) = signedInertia(mcI(rotorMass, largeRotorCOM, largeRotorRotInertiaY), legID); % first rotor
  model.I{idx}(13:18,13:18) = signedInertia(mcI(rotorMass, smallRotorCOM, smallRotorRotInertiaY), legID); % second rotor
  model.I{idx}(19:24,19:24) = signedInertia(mcI(ankleMass, ankleCOM, ankleRotInertia), legID); % second link (foot)
  %model.appearance.body{idx} = {{'3dobj','obj/mit_humanoid/shank.obj'},{},{},...
  %				 {'3dobj','obj/mit_humanoid/foot.obj'}};
  kcyl1 = [0 1.25*leg_rad 0;...
          0 -1.25*leg_rad 0];
  kcyl2 = [0 0 0;...
          signedVector(ankleLocation, legID)'];
  acyl1 = [0 0.5*leg_rad 0;...
           0 -0.5*leg_rad 0];
  acyl2 = [0 0 0;...
        signedVector([footToeLength; 0; -footHeight], legID)'];
  acyl3 = [0 0 0;...
        signedVector([-footHeelLength; 0; -footHeight], legID)'];
  model.appearance.body{idx} = ...
        {{'colour',[203, 121, 247]./255,...
        'cyl', kcyl1, 1.5*leg_rad,...
        'cyl', kcyl2, leg_rad}, {}, {},...
        {'colour',[18, 6, 186]./255,...
        'cyl', acyl1, 0.8*leg_rad,...
        'cyl', acyl2, 0.5*leg_rad,...
        'cyl', acyl3, 0.5*leg_rad}};
end

for armID = ["right","left"]

  % shoulderRy 
  idx = idx + 1;
  model.joint{idx} = revoluteJointWithRotor();
  model.joint{idx}.gearRatio = shoulderRyGearRatio;
  model.joint{idx}.rotorAxis = 'y';
  model.joint{idx}.jointAxis = 'y';
  model.Xtree{idx}(1:6,:) = plux(I3, signedVector(shoulderRyLocation, armID));
  model.Xtree{idx}(7:12,:) = plux(I3, signedVector(shoulderRyRotorLocation, armID));
  model.I{idx}(1:6,1:6) = signedInertia(mcI(shoulderRyMass, shoulderRyCOM, shoulderRyRotInertia), armID);
  model.I{idx}(7:12,7:12) = signedInertia(mcI(rotorMass, smallRotorCOM, smallRotorRotInertiaY), armID);
  %model.appearance.body{idx} = {{},{}};
  cyl1 = [0 3*leg_rad 0;...
        0 -3*leg_rad 0];
  cyl2 = [0 0 0;...
        signedVector(shoulderRxLocation, armID)'];
  model.appearance.body{idx} = ...
        {{'colour',[225, 225, 227]./255,...
        'cyl', cyl1, 1.35*leg_rad,...
        'cyl', cyl2, leg_rad}, {}};

  % shoulderRx
  idx = idx + 1;
  model.joint{idx} = revoluteJointWithRotor();
  model.joint{idx}.gearRatio = shoulderRxGearRatio;
  model.joint{idx}.rotorAxis = 'x';
  model.joint{idx}.jointAxis = 'x';
  model.Xtree{idx}(1:6,:) = plux(I3, signedVector(shoulderRxLocation, armID));
  model.Xtree{idx}(7:12,:) = plux(I3, signedVector(shoulderRxRotorLocation, armID));
  model.I{idx}(1:6,1:6) = signedInertia(mcI(shoulderRxMass, shoulderRxCOM, shoulderRxRotInertia), armID);
  model.I{idx}(7:12,7:12) = signedInertia(mcI(rotorMass, smallRotorCOM, smallRotorRotInertiaX), armID);
  %model.appearance.body{idx} = {{'3dobj','obj/mit_humanoid/shoulderRx.obj'},{}};
  cyl1 = [0.75*leg_rad 0 0;...
        -0.75*leg_rad 0 0];
  cyl2 = [0 0 0;...
        signedVector(shoulderRzLocation, armID)'];
  model.appearance.body{idx} = ...
        {{'colour',[225, 225, 227]./255,...
        'cyl', cyl1, 2.5*leg_rad,...
        'cyl', cyl2, leg_rad}, {}};

  % shoulderRz
  idx = idx + 1;
  model.joint{idx} = revoluteJointWithRotor();
  model.joint{idx}.gearRatio = shoulderRzGearRatio;
  model.joint{idx}.rotorAxis = 'z';
  model.joint{idx}.jointAxis = 'z';
  model.Xtree{idx}(1:6,:) = plux(I3, signedVector(shoulderRzLocation, armID));
  model.Xtree{idx}(7:12,:) = plux(I3, signedVector(shoulderRzRotorLocation, armID));
  model.I{idx}(1:6,1:6) = signedInertia(mcI(shoulderRzMass, shoulderRzCOM, shoulderRzRotInertia), armID);
  model.I{idx}(7:12,7:12) = signedInertia(mcI(rotorMass, smallRotorCOM, smallRotorRotInertiaZ), armID);
  %model.appearance.body{idx} = {{'3dobj','obj/mit_humanoid/upperArm.obj'},{}};
  cyl1 = [0 0 0.75*leg_rad;...
        0 0 -0.75*leg_rad];
  cyl2 = [0 0 0;...
        signedVector(elbowLocation, armID)'];
  model.appearance.body{idx} = ...
        {{'colour',[225, 225, 227]./255,...
        'cyl', cyl1, 1.65*leg_rad,...
        'cyl', cyl2, leg_rad}, {}};

  % elbow
  idx = idx + 1;
  model.joint{idx} = revoluteJointWithRotor();
  model.joint{idx}.gearRatio = elbowGearRatio;
  model.joint{idx}.rotorAxis = 'y';
  model.joint{idx}.jointAxis = 'y';
  model.Xtree{idx}(1:6,:) = plux(I3, signedVector(elbowLocation, armID));
  model.Xtree{idx}(7:12,:) = plux(I3, signedVector(elbowRotorLocation, armID));
  model.I{idx}(1:6,1:6) = signedInertia(mcI(elbowMass, elbowCOM, elbowRotInertia), armID);
  model.I{idx}(7:12,7:12) = signedInertia(mcI(rotorMass, smallRotorCOM, smallRotorRotInertiaY), armID);
  %model.appearance.body{idx} = {{'3dobj','obj/mit_humanoid/foreArm.obj'},{}};
  cyl1 = [0 0.75*leg_rad 0;...
        0 -0.75*leg_rad 0];
  cyl2 = [0 0 0;...
        signedVector([0; 0; -lowerArmLength], armID)'];
  model.appearance.body{idx} = ...
        {{'colour',[225, 225, 227]./255,...
        'cyl', cyl1, 1.15*leg_rad,...
        'cyl', cyl2, leg_rad}, {}};

end

%model.camera.direction = [1 0 0];
%model.camera.up = [0 0 1];

model = model.postProcessModel();

  function m = vectorToSkewMat(v)
    m = [0, -v(3), v(2);...
         v(3), 0, -v(1);...
         -v(2), v(1), 0];
  end

  function vec = matToSkewVec(m)
    vec = 0.5 * [m(3,2) - m(2,3), m(1,3) - m(3,1), m(2,1) - m(1,2)];
  end

  function P = getPseudoInertia(inertia)
    h = matToSkewVec(inertia(1:3,4:6));
    Ibar = inertia(1:3,1:3);
    m = inertia(6,6);   
    P = zeros(4,4);
    P(1:3,1:3) = 0.5 * trace(Ibar) * eye(3) - Ibar;
    P(1:3,4) = h;
    P(4,1:3) = h';
    P(4,4) = m;
  end

  function signedv = signedVector(v, side)
    assert(size(v,1)==3 && size(v,2)==1);
    switch side
      case 'right'
        signedv = v;
      case 'left'
        signedv = [v(1); -v(2); v(3)];
    end
  end

  function signedSpatialInertia = signedInertia(spatialInertia, side)
    assert(size(spatialInertia,1)==6 && size(spatialInertia,2)==6);
    switch side
      case 'right'
        % flip spatialInertia along Y axis
        % 1) get pseudo inertia
        P = getPseudoInertia(spatialInertia);
        X = eye(4);
        X(2,2) = -1;
        P = X * P * X;
        % construct signedSpatialInertia from pseudo inertia
        I = zeros(6,6);
        m = P(4,4);
        h = P(1:3, 4);
        E = P(1:3, 1:3);
        Ibar = trace(E) * eye(3) - E;
        I(1:3,1:3) = Ibar;
        I(1:3,4:6) = vectorToSkewMat(h);
        I(4:6,1:3) = vectorToSkewMat(h)';
        I(4:6,4:6) = m * eye(3);
        signedSpatialInertia = I;
      case 'left'
        signedSpatialInertia = spatialInertia;
    end
  end
end
