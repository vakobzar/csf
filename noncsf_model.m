%This code plots the single compartment 
%perfusion model, i.e., disregarding the CSF, and the grey matter 
%concentration measurements.

filename = 'TACdata.xlsx'; %Excel spreadsheet name
sheet = 1;                 %relevant sheet number in the spreadsheet
subject_label = 'subj';    %header of the subject/image label column
data_labels= {'start','end','grey', 'AIF'};  %headers of the data columns

%parameters for the cleaned-up data spreadsheet
%filename = 'TAC_matlab.xlsx';
%data_labels= {'start', 'end','GM','AIF'}; 

subject = 5318%5385;

k1=.3; % the input = output flow = k ml of blood per minute
k2=.2
V = 1.5; % Volume of tissue in 1000 ml
x0= [k1/V, k2/V]; % initialization of the minimization 

%load the data from the spreadsheet
[data, subjects] = TACfromXls (filename, sheet, subject_label, data_labels);

%extract the index of the subject and his or her data 
subject_index = find(subjects==subject);
subject_data = data{subject_index}(:,:);

%minimize the residual
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
[x,fval,exitflag,output] = fminunc(@(x)residual(x, subject_data),x0, options);
%plot the data
B= make_brain(x(1), x(2), subject_data);

figure; % new figure

p1=plot(cell2mat(subject_data(:,2)), cell2mat(subject_data(:,3)),...
     cell2mat(subject_data(:,2)), B);
p1(2).Marker = '.';
title({['Single compartment model, i.e., disregarding the CSF,',...
       ' vs grey matter measurements'], ['for subject ' num2str(subject),...
        ' (k_1/V=' num2str(x(1)), ' k2/V=' num2str(x(2)) ')' ]});
ylabel('B(t)');
xlabel('t');
legend('Measurements', 'Model');
saveas(gcf, 'noncsf_model', 'pdf')


 







 
       

