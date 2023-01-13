function  [dtau_dq, dtau_dqd] = ID_derivatives( model, q, qd, qdd )

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd) || ~iscell(qdd)
    [q, qd, qdd] = confVecToCell(model,q,qd,qdd);
end

a_grav = model.getGravity();
IC = model.I;
I = model.I;


v = {};
a = {};



for i = 1:model.NB
    
  vp = model.getParentVariable(i, v);
  ap = model.getParentVariable(i, a,-a_grav);
  
  
  [Xup{i}, S{i}, Sd{i}, v{i}, derivs{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp);

  vp = Xup{i}*vp;
  ap = Xup{i}*ap;
  
  a{i}       = ap + S{i}*qdd{i} + Sd{i}*qd{i};
  
  S_q = derivs{i}.S_q; % Tensor giving derivative of S w.r.t. q
 
  alpha = contract( S_q, qd{i});
  beta  = contract( S_q, qdd{i});
  
  Psid{i}    = crm(vp)*S{i} + alpha;
  new_part = derivs{i}.Sdotqd_q + beta + crm(v{i})*alpha;
  
  Psidd{i}    = crm(ap)*S{i} + crm(vp)*Psid{i} + new_part; 
  Upsilond{i} = Sd{i}  + Psid{i};
  
  IC{i} = I{i};
  
  BC{i} = 2*factorFunctions(IC{i},v{i});
  h{i}  = IC{i}*v{i};
  f{i}  =  IC{i}*a{i} + crf(v{i})*h{i};
end


dtau_dq  = q{1}(1)*0 + zeros(model.NV,model.NV);
dtau_dqd = q{1}(1)*0 + zeros(model.NV,model.NV);


dtau_dq  = repmat(0*q{1}(1), model.NV,model.NV);
dtau_dqd  = repmat(0*q{1}(1), model.NV,model.NV);

for i = model.NB:-1:1
  ii = model.vinds{i};
  
  tmp1= IC{i}*S{i}; 
  tmp2 = BC{i}*S{i}    + IC{i}*Upsilond{i}; 
  tmp3 = BC{i}*Psid{i} + IC{i}*Psidd{i}     + icrf(f{i})*S{i}; 
  tmp4 = BC{i}.'*S{i}; 
  
  j = i;
  while j > 0
     jj = model.vinds{j};
     dtau_dq(ii,jj)  = tmp1.'*Psidd{j}+tmp4.'*Psid{j};
     dtau_dq(jj,ii)  = S{j}.'*tmp3;

     dtau_dqd(jj,ii) = S{j}.'*tmp2;
     dtau_dqd(ii,jj) = tmp1.'*Upsilond{j}+tmp4.'*S{j};
     
     if model.parent(j) > 0
        tmp1 = Xup{j}'*tmp1;
        tmp2 = Xup{j}'*tmp2;
        tmp3 = Xup{j}'*tmp3;
        tmp4 = Xup{j}'*tmp4;
     end
     j = model.parent(j);
  end
  
  S_q = derivs{i}.S_q;
  dtau_dq(ii,ii) = dtau_dq(ii,ii) + contractT( S_q, f{i})';
  
  if model.parent(i) > 0
     p = model.parent(i);
     IC{p} = IC{p} + Xup{i}'*IC{i}*Xup{i};
     BC{p} = BC{p} + Xup{i}'*BC{i}*Xup{i};
     f{p}  = f{p}  + Xup{i}'*f{i};
  end
end

end


function M = myCell2Mat(model,S)
%     import casadi.*
%     M = SX.zeros(6,model.NV);
%     for i = 1:length(S)
%         M(:,model.vinds{i}) = S{i};
%     end
M = cell2mat(S);
end

function out = contract(S_q, vec)
    out = zeros( size(S_q,1), size(S_q,3) );
    
    for i = 1:size(S_q,3)
       out(:,i) =  S_q(:,:,i)*vec;
    end
end

function out = contractT(S_q, vec)
    out = zeros( size(S_q,2), size(S_q,3) );
    
    for i = 1:size(S_q,3)
       out(:,i) =  S_q(:,:,i)'*vec;
    end
end


