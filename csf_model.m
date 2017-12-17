% This code generates plots for a one compartment model and a two compartment
% model using arterial blood, grey matter and CSF
% concentration measurements.

filename = 'TACdata.xlsx'; % Excel spreadsheet name
sheet = 1;                 % relevant sheet number in the spreadsheet
subject_label = 'subj';    % header of the subject/image label column
%headers of the data columns
data_labels= {'start','end','CSF', 'grey', 'AIF'}; 
diseaseT = 'AD'; 
diseaseF = 'NL'; 

%parameters for the cleaned-up data spreadsheet
%filename = 'TAC_matlab.xlsx';
%data_labels= {'start', 'end','GM','AIF'}; 

%load the data from the spreadsheet
[data, subjects, disease] = TACfromXls (filename, sheet, subject_label,...
                                diseaseT, diseaseF,data_labels);
  
% initialization of the noncsf minimization 
y0(1)=.2;  % .8 flow from blood to brain k (ml of blood per minute)
y0(2)=.2 ;  % .4 flow from brain to blood 
    
% initialization of the one compartment csf minimization

x0(1)= 0.2;  % .8 flow from blood to brain k (ml of blood per minute)
x0(2)=0.4;  % .4 flow from brain to blood 
x0(3)= 0.2;  % .15 flow from CSF to brain 
    
f=figure;
diseaseMat = cell2mat(disease);
for subject = subjects'
    %extract the index of the subject and his or her data 
    subject_index = find(subjects==subject);
    subject_data = data{subject_index}(:,:);
    
    % determine the weights
    frames = cell2mat(subject_data(:,2))-cell2mat(subject_data(:,1));
    weights = frames / sum(frames);

    options = optimoptions(@fminunc,'Algorithm','quasi-newton');
    
    %NON-CSF MODEL
    if diseaseMat(subject_index)== 1
        diseaseflag = 'Alzheimers';
    else 
        diseaseflag = 'healthy';
    end
    
     subjectStr = sprintf(['%s subject %u'], diseaseflag, subject);
    
    [csf_n(subject_index, :),fval,exitflag,output] = ...
                    fminunc(@(x)residual_single(x, subject_data,weights,...
                    @make_brain_noncsf),y0, options);
    csf_noncsf= make_brain_noncsf(csf_n(subject_index, :), subject_data);
    T=cell2mat(subject_data(:,2));
     
    %plot single compartment noncsf model
 
    p1=plot(T, cell2mat(subject_data(:,5)),...
        T, cell2mat(subject_data(:,3)),...
        T, cell2mat(subject_data(:,4)),...
        T, csf_noncsf);
    p1(3).Marker = '.';

    ylabel('Concentration');
    xlabel('t');
    legend( ['AIF measurements'],...
        ['CSF measurements'],...
        ['Gray matter measurements'],... 
        ['Grey matter concentration (model)' 10 '(k_1=' num2str(csf_n(subject_index, 1)),...
                    ', k_2=' num2str(csf_n(subject_index, 2)),')']);
     title({['Single compartment non-CSF model for '], subjectStr, ...
        ['(initialization k_1=' num2str(y0(1)),...
                        ', k_2=' num2str(y0(2)), ')']});

    fname = sprintf(['/plots/%ucomp1noncsf.pdf'], subject );
    tightfig(f);
    saveas(f, [pwd, fname], 'pdf');
    
    
    

    %ONE COMPARTMENT CSF MODEL
    [csf_s(subject_index, :),fval,exitflag,output] = ...
                    fminunc(@(x)residual_single(x, subject_data,weights,...
                    @make_brain_single),x0, options);
    csf_single= make_brain_single(csf_s(subject_index, :), subject_data); 
     %plot single compartment CSF model
    p2=plot(T, cell2mat(subject_data(:,5)),...
        T, cell2mat(subject_data(:,3)),...
        T, cell2mat(subject_data(:,4)),...
        T, csf_single);
    p2(3).Marker = '.';
    title({['Single compartment CSF model for'], subjectStr,...
        ['(initialization at k_1=' num2str(x0(1)),...
                        ', k_2+k_4=' num2str(x0(2)), ... 
                        ', k_3=' num2str(x0(3)), ')' ]});
    ylabel('Concentration');
    xlabel('t');
    legend( ['AIF measurements'],...
        ['CSF measurements'],...
        ['Gray matter measurements'],... 
        ['Grey matter concentration (model)' 10 '(k_1=' num2str(csf_s(subject_index, 1)),...
                    ', k_2+k_4=' num2str(csf_s(subject_index, 2)),...
                    ', k_3=' num2str(csf_s(subject_index, 3)), ')']);
                
    fname = sprintf(['/plots/%ucomp1.pdf'], subject );
    tightfig(f);
    saveas(f, [pwd, fname], 'pdf');
    
    
    %TWO COMPARTMENT MODEL
    % initialization of the minimization 
    V = 1;    % volume of brain tissue (in liters)
    V_CSF = 1;% volume of CSF 

    k0(1)=csf_s(subject_index,1)*V;  % .8 flow from blood to brain k (ml of blood per minute)
    k0(2)=csf_s(subject_index,2)*V/2;  % .4 flow from brain to blood 
    k0(3)=0;  % .15 flow from CSF to brain 
    k0(4)=csf_s(subject_index,3)*V/4;  % .04 flow from brain to CSF (clearance)
    k0(5)=0.1;  % .06  flow from blood to CSF
    k0(6)=csf_s(subject_index,1)*V;  % .04 flow from CSF to blood

    %minimize the csf residual
    [csf_k(subject_index, :),fval,exitflag,output] = ...
                    fminunc(@(x)residual(x, subject_data,weights,...
                    @make_brain, V, V_CSF),k0, options);

    csf= make_brain(csf_k(subject_index, :), subject_data, V, V_CSF) ;
    csfB=cell2mat(csf);
    grey_m=csfB(1,:);
    csf_m=csfB(2,:);            

    
    
    
    %plot 2 compartment model
    p3=plot(T, cell2mat(subject_data(:,5)),...
        T, cell2mat(subject_data(:,3)),...
        T, cell2mat(subject_data(:,4)),...
        T, grey_m,...
        T, csf_m);
    p3(3).Marker = '.';
    title({['Two compartment CSF model for'], subjectStr,...
        ['(initialization k_1=' num2str(k0(1)),...
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
        ['CSF concentration (model)' 10 '(k_1=' num2str(csf_k(subject_index, 1)),...
                    ', k_2=' num2str(csf_k(subject_index, 2)),...
                    ', k_3=' num2str(csf_k(subject_index, 3)),...
                    ', k_4=' num2str(csf_k(subject_index, 4)),...
                    ', k_5=' num2str(csf_k(subject_index, 5)),...
                    ', k_6=' num2str(csf_k(subject_index, 6)), ')']);
    fname = sprintf(['/plots/%ucomp2.pdf'], subject );
    tightfig(f);
    saveas(f, [pwd, fname], 'pdf')
end
save('csf_model.mat','csf_n','csf_s', 'csf_k', 'disease');
 



