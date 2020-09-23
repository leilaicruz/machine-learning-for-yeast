clear all
close all
clc

% this code is for looking for interactions in Excel databases

file_network = 'Interactors DIA.xlsx';       % choose which network you want to use
range_n = 'A:A';                          % where to find the names
tab_n = 2;
[~,gene_names,~] = xlsread(file_network,tab_n,range_n);

file_search = 'STRING phys-gen distiction.xlsx';
range_s = 'A:A';
tab_s = 2;
[~,gene_search] = xlsread(file_search,tab_s,range_s);

[size_s,~] = size(gene_search);
%%
for ii=size_s:-1:1
    if isequal(gene_search(ii,1),'5270.UM01141P0')
    else
        gene_search(ii)=[];
    end
end