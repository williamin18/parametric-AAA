function [x,training_err,test_err,epoch] = TT_ALS(A,b,x,rank,tol,max_epoches,A_test,b_test,lambda)
%x: unknown vector in TT-format
%A: left hand side matrix, rows in rank-1 format
%b: right hand side vector
%lambda: regularization parameter, it is equal to sqrt(lambda) in the paper

d = length(A);
[n_samples,~] = size(A{1});
[~,m,r] = TTsizes(x);


x = TTorthogonalizeRL(x); %orthogonalize

i = 1; %index for TT-core updated at current iteration
dir = 1;%sweep direction
%stores product of cores of Ax with index larger/smaller than the current core index
[~,yr] = Ax_right(A,x,i);      
yl = cell(d,1);
yl{1} = ones(n_samples,1);

break_counter = 0;
break_limit = 5;
err_old = 100;

for epoch = 1:max_epoches

   
        %compute the product of A with the fixed TT-cores to compute the
        %partial derivative to update the current TT-core
        residual = b - multi_r1_times_TT(A,x);
        x{i} = x{i} + TTcore_Newton(yl{i},yr{i},A{i},x{i},residual,lambda);
 

        %orthogonalize to update the next TT-core and project gradient to
        %the next TT-core by left-to-right and right-to-left sweeps
        if dir
            %alternating sweep direction after reaching the right most
            %core
            if i == d-1
                dir = 0;
            end
            

            [x{i}, R_k] = qr(x{i},0);
            x{i + 1} = h2v(R_k * v2h(x{i + 1}, m(i + 1)), m(i + 1));


            xi = reshape(x{i},[r(i), m(i), r(i+1)]);
            xi = reshape(permute(xi, [2 1 3]),m(i),[]);
            Axi = A{i}*xi;
            Axi = reshape(Axi,n_samples,r(i),r(i+1));

            yl{i+1} = zeros(n_samples,r(i+1));
            for j = 1:r(i+1)
                yl{i+1}(:,j)  = sum(yl{i}.*Axi(:,:,j),2);
            end
            
            i = i + 1;

        else
            if i == 2
                dir = 1;
            end
            
            [Q_k, R_k] = qr(v2h(x{i}, m(i))', 0);
            x{i} = h2v(Q_k', m(i));
            x{i-1} = x{i-1} * R_k';


            xi = reshape(x{i},[r(i), m(i), r(i+1)]);
            xi = reshape(permute(xi, [2 3 1]),m(i),[]);
            Axi = A{i}*xi;
            Axi = reshape(Axi,n_samples,r(i+1),r(i));
            yr{i-1} = zeros(n_samples,r(i));
            for j = 1:r(i)
                yr{i-1}(:,j) = sum(yr{i}.*Axi(:,:,j),2);
            end

            i = i -1;

        end
        
    
    training_err = norm(residual)/norm(b);
    r_test = multi_r1_times_TT(A_test,x) - b_test;
    % norm(r_test)
    test_err = norm(r_test)/norm(b_test);
    if test_err < tol || test_err/training_err > 4
        break
    end
    
    if   err_old-training_err < tol/1000
        break_counter = break_counter+1;
        if break_counter > break_limit
            break
        end
    else
        break_counter = 0;
    end
    err_old = training_err;
end

