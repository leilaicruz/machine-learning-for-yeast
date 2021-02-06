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

        
        