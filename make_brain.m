function x= make_brain (k, subject_data,V, V_CSF )

% This function approximates the solution of
% VC'(t) = k1 A(t) - (k2+k4)C(t) + k3 C_CSF(t)
% V_CSF C_CSF'(t)= k5 A(t) -(k3+k6) C_CSF(t)+k_4 C(t)
% given by 
% C(t) =  exp(tM) [ \int_0^t exp (-sM) ds ] v 
% by 
% x(t_i) =  exp(t_i M)[ \sum_{j=1}^i (t_j - t_{j-1}) ( A(t_j)* exp(-t_jM)+
%           + A(t_{j-1})* exp(-t_{j-1}M)) /2 ] v
% for discrete measurements, where
% x(t_i) = [C(t_i); C_CSF (t_i)],
% M = [-(k2+k4)/V, k3/V; k4/V_CSF, -(k3+k6)/V_CSF], and
% v = [k1/V; k5/V_CSF]
%
% INPUTS: k is  1x6 array containing k1 ... k6 
% and nx5 cell array representing the concentration measurements for 
% a given subject as follows 
%
% {t_{j-1}} are subject_data(:,1)
% {t_{j}} are subject_data(:,2)
% {C_CSF(t_i)} are subject_data(:,3) [These measurements are not used]
% {C(t_i)} are subject_data(:,4) [These measurements are not used]
% {A(t_j)} are subject_data(:,5)
%
% OUTPUT: 2 x n cell array containing the estimates of x(t_i) given below

n = size(subject_data,1);
subject_data(2:n,6) = subject_data(1:n-1,5);
subject_data(1,6)=subject_data(1,5);
% {A(t_{j-1})} are subject_data(:,6)

data =  num2cell(subject_data,2);
I = num2cell (1:n); 

v = [k(1)/V; k(5)/V_CSF];
M = [-(k(2)+k(4))/V, k(3)/V; k(4)/V_CSF, -(k(3)+k(6))/V_CSF];

x= cellfun(@(i)xti(M, v, data, i), I, 'UniformOutput', false);

% xti returns x(t_i), as defined above
function x = xti(M, v, data, i)

    tmp = cellfun(@(data_cell)convo_term(data_cell, data{i}{2}),...
                                        data(1:i), 'UniformOutput', false);    
    x=  -expm(data{i}{2} * M) * inv(M)* sum(cat(3,tmp{:}),3) * v;
       
%       convo_term returns the i-th term in the above sum given by 
%       (t_j - t_{j-1}) ( A(t_j)* exp(-t_j M)+ A(t_{j-1})* exp(-t_{j-1} M))
        function c = convo_term(data_cell, ti)
           
           c = (data_cell{5}*(expm(-data_cell{2}*M)-...
                                      expm(-data_cell{1}*M)));
           
        end
        
 end
end
