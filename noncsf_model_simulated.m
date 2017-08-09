%This code plots the single compartment 
%perfusion model, i.e., disregarding the CSF, and the grey matter 
%concentration measurements.

filename = 'TACdata.xlsx'; %Excel spreadsheet name
sheet = 1;                 %relevant sheet number in the spreadsheet
subject_label = 'subj';    %header of the subject/image label column
data_labels= {'start','end','CSF', 'grey', 'AIF'};  %headers of the data columns

%parameters for the cleaned-up data spreadsheet
%filename = 'TAC_matlab.xlsx';
%data_labels= {'start', 'end','GM','AIF'}; 

subject = 5318%5385;

%parameters of the simulated data
k1=.69; % the input = output flow = k ml of blood per minute
k2=.20; % Volume of tissue in 1000 ml
sigma =0.1 % noise level

%intialization of the minimzation
k1divVinit=2; % the input = output flow = k ml of blood per minute
k2divVinit = .1; % Volume of tissue in 1000 ml
x0= [k1divVinit,k2divVinit]; % initialization of the minimization 


%load the data from the spreadsheet
[data, subjects] = TACfromXls (filename, sheet, subject_label, data_labels);

%extract the index of the subject and his or her data 
subject_index = find(subjects==subject);
subject_data = data{subject_index}(:,:);

%create weights based on the lenght of the measurement intervals 
%for residual minimization
frames = cell2mat(subject_data(:,2))-cell2mat(subject_data(:,1));
weights = ones(1,size(subject_data,1));
% frames / sum(frames);

%simulate an ideal brain curve
Bideal= make_brain([k1, k2], subject_data);

%apply a pseudo-random 10% error to both blood and brain curves, 
A = cell2mat (subject_data(:,5))
[n ~] = size (A);
Apert = A + 2*sigma* A.* (rand(n,1) - .5*ones(n,1)); 
Bpert = Bideal + 2* sigma* Bideal.* (rand(1,n) - .5*ones(1,n));

subject_data(:,5) = num2cell (Apert);
subject_data(:,4) = num2cell (Bpert);


%minimize the residual
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
options.MaxFunctionEvaluations = 12000
[x,fval,exitflag,output] = fminunc(@(x)residual(x, subject_data,...
                        weights, @make_brain),x0, options);
%plot the data
B= make_brain(x, subject_data);

figure; % new figure

p1=plot(cell2mat(subject_data(:,2)), Bideal, ...
        cell2mat(subject_data(:,2)), cell2mat(subject_data(:,4)),...
        cell2mat(subject_data(:,2)), B);
p1(2).Marker = '.';
title({['Single compartment model simulations, i.e., ',...
        'disregarding the CSF,'],...
        ['for AIF from subject ' num2str(subject),...
        ' (k1_{pert}/V=' num2str(x(1)), ' k2_{pert}/V=' num2str(x(2)) ')' ]});
ylabel('B(t)');
xlabel('t');
legend(['Noiseless simulation (k1/V=' num2str(k1), ' k2/V=' num2str(k2) ')' ],...
       ['Noise perturbed (sigma=', num2str(sigma), ')'],...
       ['Recovery from Noise (k1_{init}/V=' num2str(k1divVinit),...
                            ' k2_{init}/V=' num2str(k2divVinit) ')']);
saveas(gcf, 'noncsf_model_simulated', 'pdf')


 

       

