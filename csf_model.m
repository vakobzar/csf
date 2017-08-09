% This code generates a single compartment 
% perfusion model from  arterial blood, grey matter and CSF
% concentration measurements.

filename = 'TACdata.xlsx'; % Excel spreadsheet name
sheet = 1;                 % relevant sheet number in the spreadsheet
subject_label = 'subj';    % header of the subject/image label column
%headers of the data columns
data_labels= {'start','end','CSF', 'grey', 'AIF'}; 

%parameters for the cleaned-up data spreadsheet
%filename = 'TAC_matlab.xlsx';
%data_labels= {'start', 'end','GM','AIF'}; 

subject = 5318 %5385;

k1=.2; % input flow from blood to brain k ml of blood per minute
k2=.2; % output flow from brain to blood 
k3=.2; % input flow from CSF to brain 
k4=.2; % output flow from brain to CSF (clearance)
V = 1; % Volume of tissue in 1000 ml

x0= [k1, k2, k3, k4]/V; % initialization of the minimization 

%load the data from the spreadsheet
[data, subjects] = TACfromXls (filename, sheet, subject_label, data_labels);

%extract the index of the subject and his or her data 
subject_index = find(subjects==subject);
subject_data = data{subject_index}(:,:);

frames = cell2mat(subject_data(:,2))-cell2mat(subject_data(:,1));
weights = frames / sum(frames);

%minimize the non-CSF residual
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
[noncsf_k,fval,exitflag,output] = ...
                    fminunc(@(x)residual(x, subject_data, weights,...
                    @make_brain), x0(1:2), options);
%plot the data
noncsfB= make_brain(noncsf_k, subject_data);

%minimize the csf residual
[csf_k,fval,exitflag,output] = ...
                    fminunc(@(x)residual(x, subject_data,weights,...
                    @make_brain_csf),x0, options);

%plot the data
csfB= make_brain_csf(csf_k, subject_data);

figure; % new figure
T=cell2mat(subject_data(:,2));
p1=plot(T, cell2mat(subject_data(:,5)),...
        T, cell2mat(subject_data(:,3)),...
        T, cell2mat(subject_data(:,4)),...
        T, noncsfB,...
        T, csfB);
p1(3).Marker = '.';
title({['Single compartment model CSF vs. non-CSF'], ['for subject ' num2str(subject),...
        ' (initialization k_1=' num2str(x0(1)),...
                        ' k_2=' num2str(x0(2)), ... 
                        ' k_3=' num2str(x0(3)), ... 
                        ' k_4=' num2str(x0(4)),  ')' ]});
ylabel('Concentration');
xlabel('t');
legend( ['AIF measurements'],...
        ['CSF measurements'],...
        ['Gray matter measurements'],... 
        ['Brain concentration (non-CSF model)' 10 '(k_1=' num2str(noncsf_k(1)),...
                        ' k_2=' num2str(noncsf_k(2)),')'],...
        ['Brain concentration (CSF model)' 10 '(k_1=' num2str(csf_k(1)),...
                    ' k_2=' num2str(csf_k(2)),...
                    ' k_3=' num2str(csf_k(3)),...
                    ' k_4=' num2str(csf_k(4)), ')']);
saveas(gcf, 'csf_model', 'pdf')