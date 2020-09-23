clear all
close all
clc

% this code is for looking for interactions in Excel databases

file_network = 'Overview interactions.xlsx';       % choose which network you want to use
range_n = 'B1:AQ1';                                % where to find the names
[~,gene_names,~] = xlsread(file_network,range_n);

file_search = 'Interactors RAC1 physical.xlsx';
range_s = 'C:C';
tabs = 4;

[~,size_g] = size(gene_names);
select_int = zeros(tabs,size_g);

for jj=1:tabs
    [~,all_interactors,~] = xlsread(file_search,jj,range_s);
    [size_all,~] = size(all_interactors);
    for ii=1:size_g
        xx = contains(all_interactors,gene_names(ii));
        if xx == zeros(size_all)
            select_int(jj,ii) = 0;
        else
            select_int(jj,ii) = 1;  % hier nog bijschrijven onderscheid tussen P en G
        end
    end
end