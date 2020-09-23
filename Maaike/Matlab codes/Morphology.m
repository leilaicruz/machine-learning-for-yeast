clear all
close all
clc

%% this code sorts (morphological) data, compatible with the order of gene names

file_morph = 'SCMD2 datasets.xlsx';
tab_ess = 1;
tab_non_ess = 2;
tab_wt = 3;
range = 'B:G';

ess = xlsread(file_morph,tab_ess,range);
non_ess = xlsread(file_morph,tab_non_ess,range);
wt = xlsread(file_morph,tab_wt,range);

[~,name_ess] = xlsread(file_morph,tab_ess,'A:A');
[~,name_non_ess] = xlsread(file_morph,tab_non_ess,'A:A');

mean_ess = mean(ess);
mean_non_ess = mean(non_ess);
mean_wt = mean(wt);
%%
file_names = 'Genes S. cerevisiae.xlsx';
[~,gene_names] = xlsread(file_names,'A:A');
[~,gene_codes] = xlsread(file_names,'B:B');
%%
[size_n,~] = size(gene_names);
[size_ess,~] = size(ess);
[size_non_ess,~] = size(non_ess);
size_tot = size_ess + size_non_ess;

sorted = zeros(size_n,6);

for ii=1:size_n
    for jj=1:size_ess
        if isequal(gene_codes(ii),name_ess(jj))
            sorted(ii,:) = ess(jj,:);
        end
    end
    for jj=1:size_non_ess
        if isequal(gene_codes(ii),name_non_ess(jj))
            sorted(ii,:) = non_ess(jj,:);
        end
    end
end
%%

% selecteren op essential/non-essential
% scatterplot van verschillende eigenschappen
% uiteindelijk: de aanwezigheid van dit gen heeft invloed op deze
% eigenschappen

figure

subplot(2,2,1)
semilogy([1:114],wt(:,1),'o')
hold on
semilogy([1:7],ess(:,1),'o')
semilogy([1:29],non_ess(:,1),'o')
hold off

subplot(2,2,2)
semilogy([1:114],wt(:,2),'o')
hold on
semilogy([1:7],ess(:,2),'o')
semilogy([1:29],non_ess(:,2),'o')
hold off

subplot(2,2,3)
semilogy([1:114],wt(:,3),'o')
hold on
semilogy([1:7],ess(:,3),'o')
semilogy([1:29],non_ess(:,3),'o')
hold off

subplot(2,2,4)
semilogy([1:114],wt(:,4),'o')
hold on
semilogy([1:7],ess(:,4),'o')
semilogy([1:29],non_ess(:,4),'o')
hold off