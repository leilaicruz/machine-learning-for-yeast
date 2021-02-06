% This code creates a list of all the first interactors of the core network

% INPUT
load('first_interactors') % file with all interactors of present genes from Yeastmine

% OUTPUT
% gene_names_XL_phys

% fill the empty standard name cells with systematic name
index = find(cellfun(@isempty,first_interactors));
[size_fi,~] = size(first_interactors);
first_interactors(index) = first_interactors(index-size_fi);

% create a list of physical interactors
gene_names_XL_phys = first_interactors(strcmpi(first_interactors(:,5),'physical'),:);
gene_names_XL_phys = [gene_names_XL_phys(:,4);first_interactors(:,2)];
gene_names_XL_phys = unique(gene_names_XL_phys);

% change the name of ADE5,7 as Yeastmine doesn't take it
gene_names_XL_phys(strcmp(gene_names_XL_phys,'ADE5,7')) = {'YGL234W'};

clear index
clear size_fi