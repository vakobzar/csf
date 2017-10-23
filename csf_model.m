% This code generates a two compartment 
% perfusion model from arterial blood, grey matter and CSF
% concentration measurements.

filename = 'TACdata.xlsx'; % Excel spreadsheet name
sheet = 1;                 % relevant sheet number in the spreadsheet
subject_label = 'subj';    % header of the subject/image label column
%headers of the data columns
data_labels= {'start','end','CSF', 'grey', 'AIF'}; 

%parameters for the cleaned-up data spreadsheet
%filename = 'TAC_matlab.xlsx';
%data_labels= {'start', 'end','GM','AIF'}; 

subject =  5318;
V = 1;    % volume of brain tissue (in liters)
V_CSF = 1;% volume of CSF 

% initialization of the minimization
k0(1)=.8;  % .2 flow from blood to brain k (ml of blood per minute)
k0(2)=.4;  % .1 flow from brain to blood 
k0(3)=.15;  % .1 flow from CSF to brain 
k0(4)=.04;  %  flow from brain to CSF (clearance)
k0(5)=0.3;  % .06  flow from blood to CSF
k0(6)=0.4;  % .04 flow from CSF to blood

%load the data from the spreadsheet
[data, subjects] = TACfromXls (filename, sheet, subject_label, data_labels);

%extract the index of the subject and his or her data 
subject_index = find(subjects==subject);
subject_data = data{subject_index}(:,:);

frames = cell2mat(subject_data(:,2))-cell2mat(subject_data(:,1));
weights = frames / sum(frames);

options = optimoptions(@fminunc,'Algorithm','quasi-newton');

%minimize the csf residual
[csf_k,fval,exitflag,output] = ...
                    fminunc(@(x)residual(x, subject_data,weights,...
                    @make_brain, V, V_CSF),k0, options);

%plot the data
csf= make_brain(csf_k, subject_data, V, V_CSF) ;
csfB=cell2mat(csf);
grey_m=csfB(1,:);
csf_m=csfB(2,:);

figure; % new figure
T=cell2mat(subject_data(:,2));
p1=plot(T, cell2mat(subject_data(:,5)),...
        T, cell2mat(subject_data(:,3)),...
        T, cell2mat(subject_data(:,4)),...
        T, grey_m,...
        T, csf_m);
p1(3).Marker = '.';
title({['Two compartment CSF model'], ['for subject ' num2str(subject),...
        ' (initialization k_1=' num2str(k0(1)),...
                        ', k_2=' num2str(k0(2)), ... 
                        ', k_3=' num2str(k0(3)), ... 
                        ', k_4=' num2str(k0(4)), ...
                        ', k_5=' num2str(k0(5)), ...
                        ', k_6=' num2str(k0(6)), ')' ]});
ylabel('Concentration');
xlabel('t');
legend( ['AIF measurements'],...
        ['CSF measurements'],...
        ['Gray matter measurements'],... 
        ['Grey matter concentration (model)'],...
        ['CSF concentration (model)' 10 '(k_1=' num2str(csf_k(1)),...
                    ', k_2=' num2str(csf_k(2)),...
                    ', k_3=' num2str(csf_k(3)),...
                    ', k_4=' num2str(csf_k(4)),...
                    ', k_5=' num2str(csf_k(5)),...
                    ', k_6=' num2str(csf_k(6)), ')']);
saveas(gcf, 'csf_model', 'pdf')