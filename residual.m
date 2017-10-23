function r= residual(x, subject_data, weights, fun, V, V_CSF)
% This function approximates the L2 norm distance
% between the modeled and measured concentration 
    x
    GM(:,1)= cell2mat(subject_data(:,4));
    GM(:,2)= cell2mat(subject_data(:,3));
    B = cell2mat(fun(x, subject_data, V, V_CSF));
    w(:,1)=sqrt(weights);
    w(:,2)=sqrt(weights);
 
    r = norm ((GM' - B).* w'); 
end