function C= make_brain_noncsf (k, subject_data)

% This function approximates the solution of
% C'(t) = k1 A(t) - k2C(t) 
% given by 
% C(t) =  exp(-k2t) [ \int_0^t (k1 *A(u))exp (k2 u)  du ]  
%
% C(t_i)=exp(-t_{i}*k2)*k1/k2\sum_{j=1}^i x_j 
% x_j= A(t_j) *[exp(k2*t_{j})-exp(k2*(t_{j-1})]

%
% INPUTS: k is  1x2 array containing k1 ... k2 
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

% Cti returns C(t_i), as defined above
function C = Cti(k, data, i)
    ti=data{i}{2};
    tmp = cellfun(@(data_cell)convo_term(data_cell,ti, k),data(1:i));
   
    C= exp(-ti*k(2))*k(1)/k(2)* sum(tmp) ;
       
%       convo_term returns the j-th term in the above sum given by 
%       x_j= A(t_j) *[exp(k2*t_{j})-exp(k2*(t_{j-1})]
        function c = convo_term(data_cell, ti,k)
           tjminus1 = data_cell{1};
           tj = data_cell{2};
           A=data_cell{5};
           c = A*(exp(k(2)*tj)- exp(k(2)*tjminus1));
           
        end
        
 end
end
