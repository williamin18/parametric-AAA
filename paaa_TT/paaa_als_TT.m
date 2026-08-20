function [bf,info] = paaa_als_TT(samples,sampling_values,nodes_part,coef_rank,options)
% SOLVE_ALS Solve low-rank p-AAA LS problem using ALS (Alternating Least Squares).
%
%   [BF, INFO] = SOLVE_ALS(SAMPLES, SAMPLING_VALUES, ITPL_PART, COEF_RANK, OPTIONS)
%   Solve the low-rank p-AAA LS problem by using ALS.
%
%   Inputs:
%       SAMPLES          - Tensor toolbox tensor of samples to be approximated.
%       SAMPLING_VALUES  - Cell array of sampling points in each variable.
%       NODES_PART       - Cell array of node partitions in each variable.
%       COEF_RANK        - Constraint for number of terms included in the CP decomposition used to represent the barycentric coefficients.
%       OPTIONS          - Struct containing options:
%                            * options.real_transforms - Struct with fields UL, UR for real transformations.
%                            * options.coefs_init      - Initial coefficient matrices (default: random).
%                            * options.change_tol      - Convergence tolerance for ALS based on relative change of LS errors (default: 1e-2).
%                            * options.max_iter        - Maximum ALS iterations (default: 50).
%                            * options.abs_tol         - Convergence tolerance for ALS based on objective function value (default: 0).
%                            * options.more_info       - Wether or not to save and store iterates.
%
%   Outputs:
%       BF               - Rational approximant in barycentric form.
%       INFO             - Structure containing information about iterates:
%                            * info.als_bf_iterates                 - Barycentric forms for all iterations.
%                            * info.als_rel_linearized_ls_errors    - Objective function (relative linerized LS error) for all iterations.
%                            * info.als_rel_ls_errors               - Relative nonlinear LS error for all iterations.
%                            * info.als_rel_max_errors              - Relative maximum error for all iterations.
%

orders = cellfun("length", nodes_part);
num_vars = length(sampling_values);

if ~isfield(options,'real_transforms')
    options.real_transforms = struct;
end

if ~isfield(options,'coefs_init')
    coefs = arrayfun(@(ord) randn(ord, coef_rank), orders, 'UniformOutput', false);
else
    coefs = options.coefs_init;
end

if ~isfield(options,'change_tol')
    options.change_tol = 1e-2;
end

if ~isfield(options,'max_iter')
    options.max_iter = 50;
end

if ~isfield(options,'abs_tol')
    options.abs_tol = 0;
end

if ~isfield(options,'more_info')
    options.more_info = false;
end

%init
[d,m_max,r_samples] = TTsizes(samples);
[~,m_coefs,r_coefs] = TTsizes(coefs);
%Cauchy matrix
C = cell(d,1);
for i = 1:d
    itpl_indices = nodes_part{i};
    ls_indices = setdiff(1:n, itpl_indices);
    C{i} = zeros(m_max(i),m_coefs(i));
    C{i}(itpl_indices,:) = eye(m_coefs(i));
    C{i}(ls_indices,:) = 1 ./ (sampling_values{i}(ls_part) - sampling_values{i}(itpl_part).');
    C{i} = normalize(C{i},2,"norm");
end

%error = C*(itpl_samples.*alpha) - samples.*(C*alpha). First use the initial
%guess, compute orthogonalized error TT-cores to optimize the error norm
error_TT = cell(d,1);
r_prod = r_samples.*r_coefs;
for i = 1:d-1
    samples_i = reshape(samples{i},r_samples(i),m_max(i),r_samples(i+1));
    itpl_samples_i = samples_i(:,nodes_part{i},:);
    
    %Element wise product(itpl_samples.*alpha) by Kronecker product of
    %TT-cores 
    G = reshape(itpl_samples_i,[r_samples(i),1,m_coefs(i),r_samples(i+1),1]) .* ...
        reshape(coefs{i},[1,r_coefs(i),m_coefs(i),1,r_coefs(i+1)]);
    G = reshape(G,[r_prod(i), m_coefs(i), r_prod(i+1)]);

    %Get TT-core i of C*(itpl_samples.*alpha) of size r_prod(i)*m_max(i)*r_prod(i+1)
    G = permute(tensorprod(C{i},G,2,2), [2 1 3]);

    
    %Compute (C*alpha)
    H = tensorprod(C{i} , reshape(coefs{i}, [r_coefs(i),m_coefs(i),r_coefs(i+1)]) ,2 ,2);
    H = permute(H,[2 1 3]);

    %Element wise product samples.*(C*alpha) by Kronecker product of
    %TT-cores
    H = reshape(samples_i, r_samples(i),1,m_max(i),r_samples(i+1),1) .* ...
        reshape(H,1,r_coefs(i),m_max(i),1,r_coefs(i+1));
    H = reshape(H,r_prod(i),m_max(i),r_prod(i+1));

    if i == 1
        %i = 1, r_prod(1) = 1
        error_TT{i} = [reshape(G,m_max(i),r_prod(i+1)) -reshape(H,m_max(i),r_prod(i+1))];
        [error_TT{i}, R] = qr(error_TT{i},'econ');
    else
        error_TT{i}= zeros(r_prod(i)*2,m_max(i),r_prod(i+1)*2);
        error_TT{i}(1:r_prod(i),:,1:r_prod(i+1)) = G;
        error_TT{i}(r_prod(i)+1:2*r_prod(i),:,r_prod(i+1)+1:2*r_prod(i+1)) = -H;
        error_TT{i} = R*reshape(error_TT{i},r_prod(i)*2,m_max(i)*r_prod(i+1)*2);
        [error_TT{i}, R] = qr(h2v(error_TT{i},m_max(i)));
    end
end

%Start ALS from the last TT-core
direction = -1;



% initialize difference in objective functions between iterations
diff = Inf;




L1 = als_mat(samples, sampling_values, nodes_part, coefs, 1, options.real_transforms);

norm_2_samples = norm(samples(:))^2;
prev_obj_fun = norm(L1*coefs{1}(:))^2 / norm_2_samples;
obj_fun = prev_obj_fun;

fprintf('    ALS Initial      | obj fun   %.3e\n', prev_obj_fun);

info.als_bf_iterates = {};
info.als_rel_linearized_ls_errors = [];
info.als_rel_ls_errors = [];
info.als_rel_max_errors = [];

if options.more_info
    if isfield(options.real_transforms,'UL') && isfield(options.real_transforms,'UR')
        coefs{1} = options.real_transforms.UR * coefs{1};
    end
    
    itpl_samples = samples(nodes_part{:});
    itpl_samples = tensor(itpl_samples,orders);
    
    denom_coefs = ktensor(coefs);
    num_coefs = itpl_samples .* full(denom_coefs);
    itpl_nodes = cellfun(@(sv, lp) sv(lp), sampling_values, nodes_part, 'UniformOutput', false);
    
    bf = LowRankBarycentricForm(itpl_nodes,num_coefs,denom_coefs);
    
    info.als_bf_iterates{end+1} = bf;
    info.als_rel_linearized_ls_errors(end+1) = obj_fun;
    err_mat = abs(samples-bf.eval(sampling_values));
    err_mat(isnan(err_mat)) = 0; % account for exactly zero barycentric coefficients
    info.als_rel_ls_errors = norm(err_mat(:))^2 / norm_2_samples;
    max_samples = max(abs(samples),[],'all');
    info.als_rel_max_errors = max(err_mat,[],'all') / max_samples;
end

i = 0;
while diff > options.change_tol && i <= options.max_iter && obj_fun > options.abs_tol
    i = i + 1;

    for j = 1:num_vars

        % get Loewner-type matrix for LS problem
        L = als_mat(samples, sampling_values, nodes_part, coefs, j, options.real_transforms);

        % construct second gsvd matrix for ||coef|| = 1 constraint
        ckr = ones(1,coef_rank);
        for k = 1:num_vars
            if k ~= j
                ckr = khatrirao(ckr,coefs{k});
            end
        end

        % rank-deficient constraint matrix B -> need to reduce coef_rank
        rank_ckr = rank(ckr);
        if rank_ckr < coef_rank
            fprintf('        Rank Change  | reducing coef_rank: %d → %d (GSVD rank-deficient)\n', coef_rank, rank_ckr);
            coef_rank = rank_ckr;
            coefs = cellfun(@(c)c(:,1:coef_rank),coefs,'UniformOutput',false);
            L = als_mat(samples, sampling_values, nodes_part, coefs, j, options.real_transforms);
            ckr = ckr(:,1:coef_rank);
        end

        % construct B matrix for gsvd
        B = kron(ckr,eye(orders(j)));

        % compute gsvd 
        [~,~,X,SL,SB] = gsvd(L,B,0);

        % find minimum for generalized singular values
        GS = diag(SL) ./ diag(SB);
        [~,mi_GS] = min(GS);

        % get smallest "generalized singular vector" as LS solution
        S_inv = diag(1 ./ diag(SB));

        % in rare cases X will be numerically rank deficient and can not be inverted
        coef_new = lsqminnorm(X',S_inv(:,mi_GS));

        % reshape into correct form
        coef_new = reshape(coef_new,orders(j),coef_rank);

        % unfortunately vecnorm does not always work, hence the code below
        mu = arrayfun(@(col) norm(coef_new(:, col)), 1:size(coef_new, 2));
        coef_new = coef_new ./ mu;
        coefs{j} = coef_new;
    end

    % return norm scalings my to the coefficient tensor
    coefs{1} = mu .* coefs{1};

    L1 = als_mat(samples, sampling_values, nodes_part, coefs, 1, options.real_transforms);

    obj_fun = norm(L1 * coefs{1}(:))^2 / norm_2_samples;

    diff = abs((obj_fun - prev_obj_fun) / obj_fun);

    fprintf('        ALS Iter %3d | obj fun   %.3e | change %.3e\n', i, obj_fun, diff);

    if options.more_info
        real_coefs = coefs;
        if isfield(options.real_transforms,'UL') && isfield(options.real_transforms,'UR')
            real_coefs{1} = options.real_transforms.UR * real_coefs{1};
        end
        
        itpl_samples = samples(nodes_part{:});
        itpl_samples = tensor(itpl_samples,orders);
        
        denom_coefs = ktensor(real_coefs);
        num_coefs = itpl_samples .* full(denom_coefs);
        itpl_nodes = cellfun(@(sv, lp) sv(lp), sampling_values, nodes_part, 'UniformOutput', false);
        
        bf = LowRankBarycentricForm(itpl_nodes,num_coefs,denom_coefs);
        info.als_bf_iterates{end+1} = bf;
        err_mat = abs(samples-bf.eval(sampling_values));
        info.als_rel_ls_errors(end+1) = norm(err_mat(:))^2 / norm_2_samples;
        info.als_rel_max_errors(end+1) = max(err_mat,[],'all') / max_samples;
    end

    info.als_rel_linearized_ls_errors(end+1) = obj_fun;
    
    % ALS is guaranteed to monotonically decrease the error, break if the error goes up
    if prev_obj_fun < obj_fun
        break
    end
    prev_obj_fun = obj_fun;
  
end

if isfield(options.real_transforms,'UL') && isfield(options.real_transforms,'UR')
    coefs{1} = options.real_transforms.UR * coefs{1};
end

itpl_samples = samples(nodes_part{:});
itpl_samples = tensor(itpl_samples,orders);

denom_coefs = ktensor(coefs);
num_coefs = itpl_samples .* full(denom_coefs);
itpl_nodes = cellfun(@(sv, lp) sv(lp), sampling_values, nodes_part, 'UniformOutput', false);

bf = LowRankBarycentricForm(itpl_nodes,num_coefs,denom_coefs);

end
