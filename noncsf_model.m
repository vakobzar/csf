%This program plots the difference between the single compartment 
%perfusion model, i.e., disregarding the CSF, and the grey matter 
%concentration measurements


filename = 'TACdata.xlsx'; %Excel spreadsheet name
sheet = 1;                 %relevant sheet number in the spreadsheet
subject_label = 'subj';    %header of the subject/image label column
data_labels= {'start','end','grey', 'AIF'};  %headers of the data columns

%parameters for the cleaned-up data spreadsheet
%filename = 'TAC_matlab.xlsx';
%data_labels= {'start', 'end','GM','AIF'}; 

subject = 5385;

k=.3; % the input = output flow = k ml of blood per minute
V = 1.5; % Volume of tissue in 1000 ml
x0= [k,V]; % initialization of the minimization 

%load the data from the spreadsheet
[data, subjects] = TACfromXls (filename, sheet, subject_label, data_labels);

%extract the index of the subject and his or her data 
subject_index = find(subjects==subject);
subject_data = data{subject_index}(:,:);

%minimize the residual
[x,fval,exitflag,output] = fminunc(@(x)residual(x, subject_data),x0,options);
%plot the data
B= make_brain(x(1), x(2), subject_data);

figure; % new figure

p1=plot(cell2mat(subject_data(:,2)), cell2mat(subject_data(:,3)),...
     cell2mat(subject_data(:,2)), B1);
p1(2).Marker = '.';
title({['Single compartment model, i.e., disregarding the CSF,',...
       ' vs grey matter measurements'], ['for subject ' num2str(subject),...
        ' (k=' num2str(x(1)), ' V=' num2str(x(2)) ')' ]});
ylabel('B(t)');
xlabel('t');
legend('Measurements', 'Model');
saveas(gcf, 'nonCSFmodel', 'pdf')

function r= residual (x, subject_data)
% This function approximates the L2 norm distance
% between the concentration measurements of the 
    k=x(1);
    V=x(2);
    B = make_brain(k, V, subject_data);
    r = norm (cell2mat(subject_data(:,3))- B,2); 
end
 
function B= make_brain (k, V, subject_data)

% This function approximates the solution of
% VC'(t) = kA(t) - kC(t)
% given by 
% C(t) = 1/V * \int_0^t A(u) exp (-k(t-u)/V) du 
% as 
% C(t_i)=1/V*\sum_{j=1}^i \int_{t_{j-1}}^{t_j} A(t_j) exp(-k(t-u)/V)du
%       = 1/k*\sum_{j=1}^i A(t_j) [exp(-k(t-u)/V)]_{u=t_{j-1}}^{t_j}
% for discrete measurements, where 
% {t_{j-1}} are subject_data(:,1)
% {t_{j}} are subject_data(:,2)
% {C(t_i)} are subject_data(:,3)
% {A(t_j)} are subject_data(:,4)

A =  num2cell(subject_data,2);
[n,~] = size(A);
lin_idx = 1:n;
I =  num2cell (lin_idx); 
B= cellfun(@(i)sum_convo_terms(k, V, A, i), I); 

 function s = sum_convo_terms(k,V, A, i)
       C = cellfun(@(a)convo_term(a, A{i}{2}), A(1:i));  
       s = sum(C)/k ;
       
       function c = convo_term(Aj, ti)
          c = Aj{4}*(exp(-k*(ti-Aj{2})/V)-exp (-k*(ti-Aj{1})/V));
       end
 end
end


function [data, subjects]= TACfromXls (filename, sheet, subject_label, data_labels)

% This function loads the content of sheet number 'sheet' from an Excel 
% file 'filename' into a 1 x m cell array, where m is the number of
% subjects identified by a unique unsigned integer in the column headed
% by 'subject_label'.  Each cell contains a n x k double cell array where 
% n is the number of measurements per subject and k is the number of 
% variables that were measured. 

[~,~,raw]= xlsread(filename,1);
txt = cellfun(@num2str,raw,'UniformOutput',0);

%identify the column 'label_col' with subject/image label
TF = contains(txt, subject_label, 'Ignorecase', true);
[~, columns] = size(txt(TF));
if columns==1 
    [label_row,label_column,v] = find(TF);
    fprintf(['The labels %s are in column %i ',...
              'of %s.\n'], subject_label, ...
                     label_column(1), filename);
    label_col=label_column(1); 
else 
    error ('Error. Cannot identify the column with the subject ',...
             'labels %s in %s.\n', subject_label,  filename);
end

%identify the subjects labeled by unique unsigned integers 'subjects'
index=cellfun(@(s)sscanf(s,'%u',1), txt(:,label_col), 'uniformoutput',false);
index_mat= cell2mat(index(~cellfun(@isempty, index)));
[~,idx]=unique(index_mat,'rows');
subjects = index_mat(idx,:);
[total_subjects, ~] = size(subjects); 
if (total_subjects>0)
    fprintf(['%u unique subject labels found in column %u ',...
              'of %s.\n'], total_subjects, label_col, filename);
else 
    error ('Error. Cannot identify the subject labels ',...
             'in column %u of %s.\n', subject_label,  filename);
end

%identify the columns 'data_col' with the data entries
data_col = [];
[~, data_fields] = size(data_labels); 
for i = 1:data_fields
    TF = contains(txt, data_labels{i}, 'Ignorecase', true);
    [row, columns] = size(txt(TF));
    if (columns == 1) && (row >=1)  
        [data_row, data_column,v] = find(TF);
        fprintf(['The entries %s are in column %i ',...
            'of %s.\n'], data_labels{i}, data_column(1), filename);
        data_col(i) = data_column(1);
    else  
        error (['Error. Cannot identify the columns with the entries ',...
           'labels %s in %s.\n'], data_labels{i}, filename);
    end 
end



for i=1:total_subjects
    TF = contains(txt(:, label_col), string(subjects(i)), 'Ignorecase', true);
    tmp = raw(TF, data_col); 
    numericTF= cellfun(@isnumeric,tmp);
    tmp(~numericTF)={nan};
    tmp(any(cellfun(@(x) any(isnan(x)),tmp),2),:) = [];
    [entries_loaded, ~] = size(tmp); 
    fprintf(['%u entries points loaded for subject %i',...
              'from %s.\n'], entries_loaded, subjects(i), filename);
    data{i} = tmp;

end 

end





 
       

