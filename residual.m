function r= residual (x, subject_data)
% This function approximates the L2 norm distance
% between the modeled and measured concentration 
    k1divV=x(1);
    k2divV=x(2);
    B = make_brain(k1divV, k2divV, subject_data);
    r = norm (cell2mat(subject_data(:,3))- B,2); 
end