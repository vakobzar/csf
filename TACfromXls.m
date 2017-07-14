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
