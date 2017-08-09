function C= make_brain_csf (k, subject_data)

% This function approximates the solution of
% VC'(t) = (k1A(t) - k2C(t)) + (k3C_CSF(t) - k4C(t))
% given by 
% C(t) =  \int_0^t (k1/V * A(u)+ k3/V C_CSF(t) exp (-(k2+k4)(t-u)/V) du 
% as 
% C(t_i)=1/(k2+k4)* \sum_{j=1}^i k1 A(t_j)+k3 C_CSF(t_j)
%         * [exp(-(k2+k4)(t_i-u)/V)]_{u=t_{j-1}}^{t_j}
% for discrete measurements. 
%
% INPUTS: 1x4 array representing k1/V, ..., k4/V below
% and nx5 cell array representing the concentration measurements for 
% a given subject as follows 
%
% {t_{j-1}} are subject_data(:,1)
% {t_{j}} are subject_data(:,2)
% {C_CSF(t_i)} are subject_data(:,3)
% {C(t_i)} are subject_data(:,4) [These measurements are not used]
% {A(t_j)} are subject_data(:,5)
%
% OUTPUT: 1 x n cell array containing the estimates of C(t_i) given below


data =  num2cell(subject_data,2);
n = size(data,1);
I = num2cell (1:n); 

C= cellfun(@(i)Cti(k, data, i), I); 

% C_ti returns C(t_i), as defined above
function C = Cti (k, data, i)
    tmp = cellfun(@(data_cell)convo_term(data_cell, data{i}{2}), data(1:i));  
    C = sum(tmp)/(k(2)+k(4));
       
%       convo_term returns the i-th term in the above sum given by 
%       k1 A(t_j) +k3 C_CSF(t_j) * [exp(-(k2+k4)(t_i-u)/V)]_{u=t_{j-1}}^{t_j}
        function c = convo_term(data_cell, ti)
            c = (k(1)*data_cell{5}+k(3)*data_cell{3})...
                 *(exp(-(k(2)+k(4)) *(ti-data_cell{2}))-...
                   exp(-(k(2)+k(4)) *(ti-data_cell{1})));
        end
 end
end
