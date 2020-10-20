% this code sorts (morphological) data, compatible with the order of gene names

file_morph = 'SCMD2 datasets.xlsx';
tab_ess = 1;
tab_non_ess = 2;
tab_wt_e = 3;
tab_wt_ne = 4;
range = 'C:H';

ess = xlsread(file_morph,tab_ess,range);
non_ess = xlsread(file_morph,tab_non_ess,range);
wt_e = xlsread(file_morph,tab_wt_e,range);
wt_ne = xlsread(file_morph,tab_wt_ne,range);

[~,name_ess] = xlsread(file_morph,tab_ess,'A:A');
[~,name_non_ess] = xlsread(file_morph,tab_non_ess,'A:A');

mean_ess = mean(ess);
mean_non_ess = mean(non_ess);
mean_wt_e = mean(wt_e);
mean_wt_ne = mean(wt_ne);
%%
figure(1)
box_index = [repmat({'wt_e'},length(wt_e),1);repmat({'wt_ne'},length(wt_ne),1);repmat({'ess'},length(ess),1);repmat({'non_ess'},29,1)];
for ii=1:6
    subplot(2,3,ii)
    boxplot([wt_e(:,ii);wt_ne(:,ii);ess(:,ii);non_ess(:,ii)],box_index)
end

figure(2)
for ii=1:6
    subplot(2,3,ii)
    hold on
    scatter(ones(114,1),wt_e(:,ii))
    scatter(2*ones(122,1),wt_ne(:,ii))
    scatter(3*ones(7,1),ess(:,ii))
    scatter(4*ones(29,1),non_ess(:,ii))
    hold off
end
%%
[h_e,p_e] = ttest2(wt_e,ess);
[h_ne,p_ne] = ttest2(wt_ne,non_ess);

clear tab_ess
clear tab_non_ess
clear tab_wt_e
clear tab_wt_ne
clear file_morph
clear range
%%
file_names = 'Genes S. cerevisiae.xlsx';
[~,gene_names] = xlsread(file_names,'A:A');
[~,gene_codes] = xlsread(file_names,'B:B');

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