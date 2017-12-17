function C= make_brain_single (k, subject_data)

% This function approximates the solution of
% VC'(t) = k1 A(t) - kC(t) + k3 C_CSF(t)
% given by 
% C(t) =  exp(-kt) [ \int_0^t (k1 *A(u) +k3 *C_CSF(u))exp (k u)  du ]  
%
% C(t_i)=exp(-t_{i}*k)/k\sum_{j=1}^i x_j 
% x_j= (k1*A(t_j)+ k3*C_CSF(t_{j}))
%       *[exp(k*t_{j})-exp(k*(t_{j-1})]
% for discrete measurements where k= k2+k4
%
% INPUTS: k is  1x6 array containing k1 ... k6 
% and nx5 cell array representing the concentration measurements for 
% a given subject as follows 
%
% {t_{j-1}} are subject_data(:,1)
% {t_{j}} are subject_data(:,2)
% {C_CSF(t_i)} are subject_data(:,3)
% {C(t_i)} are subject_data(:,4) [These measurements are not used]
% {A(t_j)} are subject_data(:,5)
%
% OUTPUT: 1 x n cell array containing the estimates of x(t_i) given below


% {A(t_{j-1})} are subject_data(:,6)

data =  num2cell(subject_data,2);
n = size(data,1);
I = num2cell (1:n); 

C= cellfun(@(i)Cti(k, data, i), I);

% Cti returns C(t_i), as defined above
function C = Cti(k, data, i)
    ti=data{i}{2};
    tmp = cellfun(@(data_cell)convo_term(data_cell,ti, k),data(1:i));
   
    C= exp(-ti*k(2))/k(2)* sum(tmp) ;
       
%       convo_term returns the j-th term in the above sum given by 
%       % x_j= (k1*A(t_j)+ k3*C_CSF(t_{j}))
%       *[exp(k/V*(t_{j})-exp(k/V*(t_{j-1}))]
        function c = convo_term(data_cell, ti,k)
           tjminus1 = data_cell{1};
           tj = data_cell{2};
           C_CSF = data_cell{3};
           A=data_cell{5};
           c = (k(1)*A+k(3)*C_CSF)*(exp(k(2)*tj)- exp(k(2)*tjminus1));
           
        end
        
 end
end
