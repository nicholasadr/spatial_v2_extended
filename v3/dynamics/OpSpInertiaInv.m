function Lambda_inv = OpSpInertiaInv(model,q, op_sp_xforms)
% Implementation of the extended force propagator algorithm (EFPA)
% Wensing, Featherstone, and Orin ICRA 2012
    

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd)
    [q] = confVecToCell(model,q);
end

%% First forward pass -- kinematics
for i = 1:model.NB
  [Xup{i}, S{i}] = model.joint{i}.kinematics(model.Xtree{i},q{i});
  IA{i} = model.I{i};
  ee_support{i} = [];
  gca{i} = [];
end

%% Book keeping setup

ee_inds = {}; % row/col indices in Lambda_inv for each end effector
ee_ind  = 0;  % counter
ee_list = []; % List of bodies with end effectors
ee_dim  = {}; % Dimension of end-effector op-sp. (e.g., point -> dim=3 ) 

for k1 = 1:length(op_sp_xforms)
    Xend = op_sp_xforms{k1};
    
    if size(Xend,1) == 0
        continue
    end
    ee_list(end+1) = k1;
    ee_dim{k1} = size(Xend,1);
    
    ee_inds{k1} = (ee_ind+1) : (ee_ind + ee_dim{k1}); % rows/cols of Op-Sp inertia for this one
    ee_ind = ee_ind + ee_dim{k1};                               % increment counter
    
    % Set up supporting bodies
    supporting_bodies = model.ancestors{k1};
    for i = supporting_bodies  % for each support body      
        ee_support{i}(end+1) = k1;
    end
    
    % Extended force propagator to end effector (just a fixed transform)
    ChiUpExtended{k1,k1} = Xend;
    
    % Determine the ancestor body where we will compute off diagonal terms
    for k2 = ee_list % for each previous op sp body        
        ancestor = max( intersect( model.ancestors{k2}, supporting_bodies) );
        gca{ancestor}(end+1,:) = [k1, k2]; % Note we also have the case when ancestor=k1=k2 here [different from ICRA alg]
    end
end

%% Main Back Pass - ABA Step
for i = model.NB:-1:1
  p = model.parent(i);
  
  U{i}  = IA{i} * S{i};
  d{i}  = S{i}.' * U{i};
  U{i}  = Xup{i}.' * U{i};

  % Force propagators (a.k.a., articulated transforms)
  % Analagous to Xup
  Chiup{i} =  (Xup{i} - (S{i}/d{i})*U{i}.');
      
  K{i} = (S{i}/d{i})*S{i}.';
  
  if p ~= 0
    Ia = Xup{i}.' * IA{i} * Xup{i} - U{i}/d{i}*U{i}.';
    IA{p} = IA{p} + Ia;
    for k = ee_support{i}
        ChiUpExtended{k,p} = ChiUpExtended{k,i}*Chiup{i};
    end
  end
end

Lambda_inv_tmp = {}; % Temp, storing acceleration at bodies caused by end effector forces
Lambda_inv = zeros(ee_ind); % Final, storing end-effector accelerations cause by body forces

%% Main forward pass -- compute the operational space inertia, finally
for i = 1:model.NB
    
    for k = ee_support{i} % foreach end-effector supported by body i
        if model.parent(i) > 0
            Lambda_inv_prev = Lambda_inv_tmp{model.parent(i) , k};
        else
            Lambda_inv_prev = zeros( 6, ee_dim{k} );
        end
        Lambda_inv_tmp{i,k} = Chiup{i}*Lambda_inv_prev ... % part due to previous body accelerating
                                  + K{i}*ChiUpExtended{k,i}.'; % part due to joint acceleration
    end
    
    gca_set = gca{i};
    for pair_ind = 1:size( gca_set, 1)
        k1 = gca_set(pair_ind, 1);
        k2 = gca_set(pair_ind, 2);
        inds1 = ee_inds{k1};
        inds2 = ee_inds{k2};
        Lambda_inv(inds1, inds2) = ChiUpExtended{k1,i}*Lambda_inv_tmp{i,k2};
        if k1 ~= k2
            Lambda_inv(inds2, inds1) = Lambda_inv(inds1, inds2).';
        end
    end
end