% This code builds an interaction matrix from all genes, present and
% absent and their first interactors.

% INPUT
load('all_first_interactors')   % file with all first interactors from Yeastmine
load('gene_names_XL')           % retrieved from A_Find_first_interactors
load('non_SC_genes');           % names of absent genes
load('non_SC_interactions');    % file with first interactors of absent genes

% OUTPUT
% interact_matrix_XXL
% gene_names_XXL

% fill the empty standard name cells with systematic name
index_afi = find(cellfun(@isempty,all_first_interactors));
[size_afi,~] = size(all_first_interactors);
all_first_interactors(index_afi)=all_first_interactors(index_afi-size_afi);

% divide the interactions in 'physical' and 'genetic'
afi_phys_i = strcmp(all_first_interactors(:,5),'physical');
afi_phys = all_first_interactors(afi_phys_i,:);
afi_gen_i = strcmp(all_first_interactors(:,5),'genetic');
afi_gen = all_first_interactors(afi_gen_i,:);

[size_XL,~] = size(gene_names_XL);
interact_matrix_XL_phys = zeros(size_XL);
interact_matrix_XL_gen = zeros(size_XL);

% create a physical interaction matrix
for ii=1:length(gene_names_XL)
    index = strcmp(afi_phys(:,2),gene_names_XL(ii));
    list = afi_phys(index,4);
    [~,pos] = intersect(gene_names_XL,list);
    interact_matrix_XL_phys(ii,pos)=1;
    clear index
    clear list
end

% create a genetic interaction matrix
for ii=1:length(gene_names_XL)
    index = strcmp(afi_gen(:,2),gene_names_XL(ii));
    list = afi_gen(index,4);
    [~,pos] = intersect(gene_names_XL,list);
    interact_matrix_XL_gen(ii,pos)=2;
    clear index
    clear list
end

% find interactions of absent genes
afi_phys_i_abs = ~strcmp(non_SC_interactions(:,3),'genetic');
afi_phys_abs = non_SC_interactions(afi_phys_i_abs,:);
afi_gen_i_abs = ~strcmp(non_SC_interactions(:,3),'physical');
afi_gen_abs = non_SC_interactions(afi_gen_i_abs,:);

size_XXL = size_XL + length(non_SC_genes);
interact_matrix_XXL_phys = interact_matrix_XL_phys;
interact_matrix_XXL_gen = interact_matrix_XL_gen;

% add the physical interactions of absent genes to the network
for ii=1:length(non_SC_genes)
    index = strcmp(afi_phys_abs(:,1),non_SC_genes(ii));
    list = afi_phys_abs(index,2);
    [~,pos] = intersect(gene_names_XL,list);
    interact_matrix_XXL_phys(pos,size_XL+ii)=1;
    interact_matrix_XXL_phys(size_XL+ii,pos)=1;
    clear index
    clear list
end

% add the physical interactions between absent genes to the network
absent_matrix_phys = [0,1,0,1,1;1,0,1,1,1;0,1,0,0,0;1,1,0,0,0;1,1,0,0,0];
interact_matrix_XXL_phys(size_XL+1:size_XXL,size_XL+1:size_XXL) = absent_matrix_phys;

% add the genetic interactions of absent genes to the network
for ii=1:length(non_SC_genes)
    index = strcmp(afi_gen_abs(:,1),non_SC_genes(ii));
    list = afi_gen_abs(index,2);
    [~,pos] = intersect(gene_names_XL,list);
    interact_matrix_XXL_gen(size_XL+ii,pos)=2;
    interact_matrix_XXL_gen(pos,size_XL+ii)=2;
    clear index
    clear list
end

% add the genetic interactions between absent genes to the network
absent_matrix_gen = [0,2,0,0,0;2,0,2,2,2;0,2,0,0,0;0,2,0,0,0;0,2,0,0,0];
interact_matrix_XXL_gen(size_XL+1:size_XXL,size_XL+1:size_XXL) = absent_matrix_gen;

% combine all information
gene_names_XXL = [gene_names_XL;non_SC_genes];
interact_matrix_XXL = interact_matrix_XXL_phys + interact_matrix_XXL_gen;

clear afi_gen
clear afi_gen_i
clear afi_phys
clear afi_phys_i
clear index_afi
clear size_afi
clear afi_gen_abs
clear afi_gen_i_abs
clear afi_phys_abs
clear afi_phys_i_abs
clear ii
clear pos
clear size_XL
clear absent_matrix_gen
clear absent_matrix_phys
clear size_XXL
clear interact_matrix_XL_gen
clear interact_matrix_XL_phys
clear interact_matrix_XXL_gen