function r= residual(x, subject_data, weights, fun)
% This function approximates the L2 norm distance
% between the modeled and measured concentration 
   
    GM= cell2mat(subject_data(:,4));
    B =  fun(x, subject_data);
    r = norm ((GM - B).*  sqrt(weights),2); 
end