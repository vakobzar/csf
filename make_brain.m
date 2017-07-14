function B= make_brain (k1divV, k2divV, subject_data)

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
% {C(t_i)} are subject_data(:,3)
% {A(t_j)} are subject_data(:,4)

A =  num2cell(subject_data,2);
[n,~] = size(A);
lin_idx = 1:n;
I =  num2cell (lin_idx); 
B= cellfun(@(i)sum_convo_terms(k1divV, k2divV, A, i), I); 

 function s = sum_convo_terms(k1divV,k2divV, A, i)
       C = cellfun(@(a)convo_term(a, A{i}{2}), A(1:i));  
       s = k1divV*sum(C)/k2divV ;
       
       function c = convo_term(Aj, ti)
          c = Aj{4}*(exp(-k2divV*(ti-Aj{2}))-exp (-k2divV*(ti-Aj{1})));
       end
 end
end
