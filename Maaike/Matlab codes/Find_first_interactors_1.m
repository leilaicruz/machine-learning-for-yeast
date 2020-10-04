% INPUT
% file with all interactors of core genes from Yeastmine

% OUTPUT
% list of names of all first interactors

first_interactors = readtable('First Interactors.csv');
first_interactors = table2cell(first_interactors);

% fill the empty standard name cells with systematic name
index_fi = find(cellfun(@isempty,first_interactors));
[size_fi,~] = size(first_interactors);
first_interactors(index_fi)=first_interactors(index_fi-size_fi);

% create a list of names
gene_names_XL = unique(first_interactors(:,4));

% change the name of ADE5,7 as Yeastmine doesn't take it
gene_names_XL(strcmp(gene_names_XL,'ADE5,7')) = {'YGL234W'};

clear index_fi
clear size_fi