function out = fcn_import_params_solidworks(filePrefix)
fileList = dir(filePrefix);       % fileList = dir('inertial_params*.txt');
nFiles = numel(fileList);

% iterate for multiple files with same prefix
for kk = 1:nFiles


    % read each data file.
    C = readcell(fileList(kk).name,'Delimiter','\r\n');     % read txt


    % detect where mass is
    idx = startsWith(C,'Mass = ');
    mass_chr = regexp(C{idx},'[-|\d|.]+','match','once');
    mass_val = str2num(mass_chr);


    % validate mass
    if isnumeric(mass_val)
        mass = mass_val*1E-3;
    else
        error('mass value is missing')
        return;
    end


    % detect where com position is
    com_ = [];
    startPattern = {'X = ', 'Y = ', 'Z = '};
    for ii = 1:numel(startPattern)
        idx = startsWith(C,startPattern{ii});
        pos_chr = regexp(C{idx},'[-|\d|.]+','match','once');
        pos_val = str2num(pos_chr);
        com_ = [com_ ; pos_val];
    end


    % validate com position
    if numel(com_) == 3
        com = com_*1E-3;
    else
        error('com position is not 3-dim')
        return;
    end


    % detect where inertia values are located
    iner_ = [];
    startPattern = {'Ixx = ', 'Iyx = ', 'Izx = '};

    for ii = 1:numel(startPattern)
        idx = startsWith(C,startPattern{ii});
        a = regexp(C{idx},'[-|\d|.]+','match');
        a_ = cell2mat(cellfun(@str2num, a, 'UniformOutput', false));

        iner_ = [iner_;a_]; 
    end


    % validate by symmetry
    if issymmetric(iner_)
        iner = iner_*1E-9;
    else
        error('failed to check inertia mtx symmetry')
        return;
    end


    % store outputs
    out.name{kk} = fileList(kk).name;
    out.mass{kk}=  mass;
    out.CoM{kk} = com;
    out.Ic{kk} = iner;

end

end
