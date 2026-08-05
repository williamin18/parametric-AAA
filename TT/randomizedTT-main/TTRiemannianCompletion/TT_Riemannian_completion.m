function [x,training_err,test_err,epoch] = TT_Riemannian_completion(sample_indices,x_samples,x,rank,tol,max_epoches,A_test,b_test)
%TT_NEWTON_GD Summary of this function goes here
%   Detailed explanation goes here
A = sample_indices;

d = length(A);
[n_samples,~] = size(A{1});
[~,m,~] = TTsizes(x);

x = TTorthogonalizeLR(x);

r = x_samples - multi_r1_times_TT(A,x);
dx_TT = 0;
beta = 0;

break_counter = 0;
break_limit = 5;
err_old = 100;

for epoch = 1:max_epoches


    [V,dUx,dx_TT] = TT_Riemannian_completion_Gradient(A,x,r,beta,dx_TT);
    x = TT_Riemannian_update(x,V,dUx,1,rank);
    
    

    r = x_samples - multi_r1_times_TT(A,x);
    beta = 1;
    
    

    training_err = norm(r)/norm(x_samples);
    r_test = multi_r1_times_TT(A_test,x) - b_test;
    test_err = norm(r_test)/norm(b_test);

    if test_err < tol
        break
    end

    if   err_old-training_err < tol/10
        break_counter = break_counter+1;
        if break_counter > break_limit
            break
        end
    else
        break_counter = 0;
    end
    err_old = training_err;
end
end

