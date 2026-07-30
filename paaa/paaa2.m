function [bf,info] = paaa2(samples,sampling_values,tol,options)

if nargin < 4
    options = struct;
end
if ~isfield(options,'max_iter')
    options.max_iter = length(samples);
end

[n_samples, num_vars] = length(sampling_values);


err_mat = abs(samples-mean(samples,'all'));
max_samples = max(abs(samples),[],'all');
[max_err,max_idx] = max(err_mat,[],'all');

j = 0;
nodes = zeros(n_samples,num_vars);
nodes_sample = zeros(n_samples);
while (max_err > max_samples * tol && j < options.max_iter)
    nodes(j,:) = sampling_values(max_idx,:);
    nodes_sample(j) = samples(max_idx,:);

    nodes_j = nodes(1:j,:);
    nodes_sample_j = nodes_sample(1:j);

    sampling_values(max_idx,:) = [];
    samples(max_idx) = [];

    C = sampling_values(:,1) - nodes_j(:,1).';
    for i = 2:num_vars
        C = C.*(sampling_values(:,i) - nodes_j(:,i).');
    end
    C = 1./C;
    L = samples.*C - (nodes_sample.').*C;

    [~,~,X] = svd(L,0);
    alpha = X(:,end);

    samples_predict = (C*(alpha.*nodes_sample))./(C*alpha);
    err_mat = abs(samples-samples_predict);
    [max_err,max_idx] = max(err_mat,[],'all');

    fprintf('p-AAA    | rel max err %.3e | rel LS err %.3e\n', max_err/max_samples, norm(err_mat)/norm(samples));

end

end

