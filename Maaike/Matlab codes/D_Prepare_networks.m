% This code prepares networks for analysis. It deletes genetic
% interactions, removes absent genes, removes self-interactions and removes
% genes which have no interactions at all.

% INPUT
load('core_names');                     % retrieved from article
load('interact_matrix_small_phys');     % retrieved from B_Build_interaction_matrix
load('existing');                       % choose existing(from Repertoire_edited) or non_existing(_0) (from C_Build_non_existing_networks)

% OUTPUT
% network
% names

ex = non_existing;                      % (non_)existing

network = cell(1,size(ex,1));
names = cell(1,size(ex,1));

for ii=1:size(ex,1)
    
    % import files with network, names and presence
    network{ii} = interact_matrix_small;
    names{ii} = core_names;

    % only look at the physical interactions    
    index = network{ii}==2;
    network{ii}(index) = 0;
    clear index

    % delete all interactions with absent genes
    index = ex(ii,:)==0;
    delete = zeros(sum(index),1);
    
    for jj=1:sum(index)
        temp = core_names(index);
        delete(jj) = find(strcmp(names{ii},temp(jj)));
    end
    
    network{ii}(delete,:)=[];
    network{ii}(:,delete)=[];
    names{ii}(delete)=[];

    clear index
    clear temp
    clear delete

    % make matrix of zeros and ones
    network{ii} = network{ii}~=0;

    % delete self interactions
    for jj=1:size(network{ii})
        network{ii}(jj,jj)=0;
    end

    % delete genes which have no interactions
    index = find(sum(network{ii})==0);
    network{ii}(index,:) = [];
    network{ii}(:,index) = [];
    names{ii}(index) = [];
    clear index

end

clear ii
clear jj
clear ex