% This code explores characteristics of individual nodes.

% INPUT
load('network_graph')
load('names')
load('core_names')
load('GTPase')
load('GAP')
load('GEF')

% OUTPUT
% degree_centrality
% closeness_centrality
% betweenness_centrality
% eccentricity_centrality
% CC

degree_centrality = cell(size(network_graph));
closeness_centrality = cell(size(network_graph));
betweenness_centrality = cell(size(network_graph));
eccentricity_centrality = cell(size(network_graph));
CC = cell(size(network_graph));
all_cent = cell(size(network_graph));
cent_names = ["betweenness","closeness","degree","eccentricity","clustering"];

for ii=1:length(network_graph)
    
    % total number of shortest paths that go through N / total number of shortest paths
    betweenness_centrality{ii} = centrality(network_graph{ii},'betweenness');
    
    % 1 / sum of distance between N and all other nodes
    closeness_centrality{ii} = centrality(network_graph{ii},'closeness');
    
    % degree of N
    degree_centrality{ii} = centrality(network_graph{ii},'degree');
        
    % 1 / max distance to another node
    dd = distances(network_graph{ii});
    dd(dd == Inf) = 0;
    eccentricity_centrality{ii} = transpose(1./max(dd));
    eccentricity_centrality{ii}(eccentricity_centrality{ii} == Inf) = 0;
    
    % 2 * edges between neighbours / degree(degree - 1)
    CC{ii} = transpose(clusteringcoef(network_graph{ii}));
    
    % gather all information in one variable
    all_cent{ii} = [betweenness_centrality{ii},closeness_centrality{ii},degree_centrality{ii},eccentricity_centrality{ii},CC{ii}];
    
end

% sort centralities per gene
count = 0;
cent_per_gene = cell(1,length(core_names));
for ii=1:length(core_names)
    for jj=1:length(all_cent)
        kk = find(strcmp(names{jj},core_names{ii}));
        if ~ isempty(kk)
            count = count+1;
            for ll=1:length(cent_names)
                cent_per_gene{ii}(count,ll) = all_cent{jj}(kk,ll);
            end
        end
    end
    count = 0;
end

% calculate averages per gene
avg_cent_e = zeros(length(cent_per_gene),length(cent_names));
for ii=1:length(cent_per_gene)
    if ~isempty(cent_per_gene{ii})
        avg_cent_e(ii,:) = sum(cent_per_gene{ii},1)./size(cent_per_gene{ii},1);
    end
end

% make distinction between GAPs, GEFs and GTPases
index_e = strings(length(core_names),1);
for ii=1:length(GTPase)
    index_e(strcmp(core_names,GTPase(ii))) = "GTPase";
end
for ii=1:length(GAP)
    index_e(strcmp(core_names,GAP(ii))) = "GAP";
end
for ii=1:length(GEF)
    index_e(strcmp(core_names,GEF(ii))) = "GEF";
end
for ii=1:length(index_e)
     if isempty(index_e{ii})
        index_e(ii) = "other";
     end
end

figure
hold on
gplotmatrix(avg_cent_e,avg_cent_e,index_e,['b' 'g' 'r' 'm'],'o',[],[],[],cent_names,cent_names);

clear jj
clear kk
clear ll
clear dd
clear count

function CC = clusteringcoef(network_graph)
% clusteringcoef    - clustering coefficient of given adjacency matrix
%
%   coefs = clusteringcoef(g) cluster coefficient is ratio of the number of
%   edges Ei between the first neighbors of the vertex i, and the
%   respective number of edges, Ei(max) = ai(ai-1)/2, in the complete graph
%   that can be formed by the nearest neighbors of this vertex:
%
%   g is a graph or an alternatively adjacency matrix.
%
%          2 Ei
%   ci = -----------
%        ai (ai - 1)
%
%   A note that do no have a link with others has clustering coefficient of NaN.


if isa(network_graph, 'graph')
    adj = adjacency(network_graph);
else
    adj = network_graph;
end

n = length(adj);
CC = zeros(1,n);
for k = 1:n
    neighbours = [find(adj(:,k))]';
    neighbours(neighbours==k) = []; % self link deleted
    a = length(neighbours);
    if a == 0; CC(k) = NaN; continue; end
    if a < 2, continue; end
    E = 0;
    for ks = 1:a
        k1 = neighbours(ks);
        for k2 = neighbours(ks+1:a)
            if adj(k1,k2) || adj(k2,k1)
                E = E + 1;
            end
        end
    end
    CC(k) = 2 * E / (a * (a-1));
end
end