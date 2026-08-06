function [y] = formRank1Tensor(ref_point,r1_samples,d)
%FORMRANK1TENSOR Summary of this function goes here
%   Detailed explanation goes here
% y = cell(d,1);
% ref_points = r1_idx(1,:);
% ref_out = r1_samples(1);
% counter = 2;
% for i = 1:d
%     temp = r1_samples(counter:(counter+m-1));
%     y{i} = [temp(1:(ref_points(i)-1)); ref_out ;temp(ref_points(i):m)];
% end
% for i = 1:d-1
%     y{i} = y{i}/ref_out;
% end
% y = TTorthogonalizeLR(y);


y = cell(d,1);
for i = 1:d-1
    y{i} = r1_samples{i}/ref_point;
end
y{d} = r1_samples{d};
y = TTorthogonalizeLR(y);
end

