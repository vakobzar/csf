function B= make_brain (k, subject_data)

% This function approximates the solution of
% VC'(t) = k1A(t) - k2C(t)
% given by 
% C(t) = k1/V * \int_0^t A(u) exp (-k2(t-u)/V) du 
% as 
% C(t_i)=k1/V*\sum_{j=1}^i \int_{t_{j-1}}^{t_j} A(t_j) exp(-k2(t-u)/V)du
%       = k1/k2*\sum_{j=1}^i A(t_j) [exp(-k2(t_i-u)/V)]_{u=t_{j-1}}^{t_j}
% for discrete measurements, where 
% {t_{j-1}} are subject_data(:,1)
% {t_{j}} are subject_data(:,2)
% {C_CSF(t_i)} are subject_data(:,3)[These measurements are not used]
% {C(t_i)} are subject_data(:,4) [These measurements are not used]
% {A(t_j)} are subject_data(:,5)

A =  num2cell(subject_data,2);
n = size(A,1);
I = num2cell (1:n); 
B=  cellfun(@(i)sum_convo_terms(k, A, i), I); 

 function s = sum_convo_terms(k, A, i)
       C = cellfun(@(a)convo_term(a, A{i}{2}), A(1:i));  
       s = k(1)*sum(C)/k(2) ;
       
       function c = convo_term(Aj, ti)
          c = Aj{5}*(exp(-k(2)*(ti-Aj{2}))-exp (-k(2)*(ti-Aj{1})));
       end
 end
end
