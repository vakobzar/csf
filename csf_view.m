hm= matfile('csf_model.mat');
csf_n=hm.csf_n;
csf_s=hm.csf_s;
csf_k=hm.csf_k;
disease=hm.disease; 

f1= figure;  
gscatter(csf_k(:, 3),csf_k(:, 4),cell2mat(disease), 'br','xo');
xlabel ('k_3');
ylabel ('k_4');
legend( ['Healthy'], ['Alzheimers']);
title ('Two compartment model with CSF');
tightfig(f1);
saveas(f1, [pwd, '/plots/twoCompClearance.pdf'], 'pdf');

f2=figure;
gscatter(csf_n(:, 2),zeros(size(csf_n,1),1), cell2mat(disease), 'br','xo');
xlabel ('k_2'); 
legend(['Healthy'],['Alzheimers'], 'Location','northeast');
ax = gca;
ax.YTick = '';
ax.YLim =[-0.03, 0.03]; 
ax.DataAspectRatio=[1 1 1];  
title ('Single compartment model without CSF');
tightfig(f2);
saveas(f2, [pwd, '/plots/SingleCompNonCSFClearance.pdf'], 'pdf');

f3=figure;
gscatter(csf_s(:, 3),zeros(size(csf_s,1),1),cell2mat(disease), 'br','xo');
xlabel ('k_2+k_4'); 
legend(['Healthy'],['Alzheimers'], 'Location','northeast');
ax = gca;
ax.YTick = '';
ax.YLim =[-0.16, 0.16]; 
ax.DataAspectRatio=[1 1 1];  
title ('Single compartment model with CSF');
tightfig(f3);
saveas(f3, [pwd, '/plots/SingleCompCSFClearance.pdf'], 'pdf');