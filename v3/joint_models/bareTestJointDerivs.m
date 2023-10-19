
%%
fprintf('=== Revolute Pair Absolute ===\n');
joint = revolutePairAbsolute();
joint.jointAxis{1} = 'x';
joint.jointAxis{2} = 'y';
joint.XtreeInternal = randomAdSE3();

q = rand(2,1);
qd = rand(2,1);
vp = rand(6,1);

Xtree = randomAdSE3();

q_derivs  = complexStepJacobian( @(x) helper(joint, Xtree, x, qd,vp), q);
qd_derivs = complexStepJacobian( @(x) helper(joint, Xtree, q, x,vp), qd);
[Xup, S, Sd, v, derivs] = joint.kinematics(Xtree, q, qd, vp);
checkEqual(derivs.Sdotqd_q, q_derivs)
checkEqual(derivs.Sdotqd_qd, qd_derivs)

%%
fprintf('=== Revolute Pair Absolute w/ Rotors===\n');

joint = revolutePairAbsoluteWithRotors();
joint.jointAxis{1} = 'x';
joint.jointAxis{2} = 'y';

joint.rotorAxis{1} = 'x';
joint.rotorAxis{2} = 'y';

joint.gearRatio{1} = 2;
joint.gearRatio{2} = 3;
joint.beltRatio{1} = 5;
joint.beltRatio{2} = 7;

joint.XtreeInternal = randomAdSE3();

q = rand(2,1);
qd = rand(2,1);
vp = rand(6,1);

Xtree = [eye(6); eye(6) ; eye(6); zeros(6)];

q_derivs  = complexStepJacobian( @(x) helper(joint, Xtree,x, qd,vp), q);
qd_derivs = complexStepJacobian( @(x) helper(joint, Xtree,q, x,vp), qd);

[Xup, S, Sd, v, derivs] = joint.kinematics(Xtree , q, qd, vp);
checkEqual(derivs.Sdotqd_q, q_derivs)
checkEqual(derivs.Sdotqd_qd, qd_derivs)



%%

function Sdot_qd = helper(joint, Xtree, q, qd, vp)
    [Xup, S, Sd, v, derivs] = joint.kinematics(Xtree, q, qd, vp);
    Sdot = Sd - crm(v)*S;
    Sdot_qd = Sdot*qd;
end

function checkEqual(a, b)
    fprintf('Equality = %e\n',norm( a(:) - b(:) ));
end

