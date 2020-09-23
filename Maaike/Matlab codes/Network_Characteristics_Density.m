clear all
close all
clc

% this code extracts which genes are present according to the Repertoire,
% matches this with Overview Interactions and fits the network to the power law
% to determine the scale-freeness and visualises it. It calculates the density, shortest paths, clustering
% coefficient, fits this to 1/k. Only apply to physical interactions.
% (builds on Network_Characteristics_Clustering_Function_Physical)

% import file of network
file_network = 'Overview interactions.xlsx';       % choose which network you want to use
original_network = xlsread(file_network,2);        % orginal_network is imported, tab 2 is without self interactions
range_n = 'B1:AQ1';                                % where to find the names
[~,gene_names,~] = xlsread(file_network,range_n);  % string of names of genes

% presence of genes
file_presence ='Repertoire_Edited.xlsx';
range_p = 'C16:AR16';                              % choose which genes are present in your network
presence = xlsread(file_presence,range_p);
[~,size_p]=size(presence);

% only look at the physical interactions
phys_network = original_network;                   
for ii=1:size_p
    for jj=1:size_p
        if phys_network(ii,jj)==2                       % legend: phys=1, gen=2, both=3
            phys_network(ii,jj)=0;
        end
    end
end

% delete all interactions with absent genes
present_network = phys_network;
for ii=size_p:-1:1                                 % backwards, otherwise for-loop stagnates
    if presence(ii)==0      
        present_network(ii,:)=[];
        present_network(:,ii)=[];
        gene_names(ii)=[];
    end
end

% make matrix of zeros and ones
binary_network = present_network./present_network;
binary_network(isnan(binary_network))=0;

% visualise the network
visualisation = graph(binary_network,gene_names);
figure
subplot(2,2,1)
plot(visualisation)

% delete genes which have no interactions
int_binary_network = binary_network;
for ii=size(int_binary_network):-1:1
    if binary_network(ii,:)==zeros(size(binary_network))
        int_binary_network(ii,:)=[];
        int_binary_network(:,ii)=[];
        gene_names(ii)=[];
    end
end

% visualise the network
int_visualisation = graph(int_binary_network,gene_names);
subplot(2,2,2)
plot(int_visualisation)

% network measures
[N_nodes,~] = size(int_binary_network);    % number of nodes is number of rows or columns
N_links = nnz(int_binary_network);         % number of links is number of non-zero elements
avg_degree = 2*N_links/N_nodes;            % average degree <k>=2L/N
degree = sum(int_binary_network);          % number of links per node

% density
N_links_max = N_nodes*(N_nodes-1)/2;
density = N_links/N_links_max;

% average pathlength calculation
pathlength = zeros(N_nodes,N_nodes);
for ii=1:N_nodes
    for jj=ii:N_nodes
        [~,length] = shortestpath(int_visualisation,ii,jj);
        pathlength(ii,jj) = length;
    end    
end
avg_pathlength = sum(pathlength,'all')/nnz(pathlength);
diameter = max(max(pathlength));

% average clustering coefficient calculation
binary_temp = int_binary_network;
nn = zeros(1,N_nodes);   % number of neighbours of neighbours
CC = zeros(1,N_nodes);   % clustering coefficient, C=2*n/k(k-1)
for jj=1:N_nodes         % determine coefficient for each node
    for ii=1:N_nodes     % turn every value for which there is no interaction zero                           
        if int_binary_network(ii,jj)==0      
            binary_temp(ii,:)=0;
            binary_temp(:,ii)=0;
        end
    end
    nn(jj) = sum(binary_temp,'all')/2;
    CC(jj)= (2*nn(jj))/(degree(jj)*(degree(jj)-1));
    CC(isnan(CC))=0;
    binary_temp=int_binary_network;
end
avg_CC = sum(CC)/N_nodes;

%degree distribution
max_degree = max(degree);
distr_count = [1:max_degree];           % nodes with no connections are not taken into account
distr_abs = histcounts(degree);         % absolute distribution
distr_rel = histcounts(degree)/N_nodes; % relative distribution

% plot of degree distribution
subplot(2,2,3)
histogram(degree,'Binwidth',5)  % which binwidth or number of bins?
xlabel('Degree')
ylabel('Counts')

% power law fit to degree distribution
[fitresult_p,gof_p] = Fit_Function_p(distr_count,distr_rel);

% 1/k fit to clustering coefficient
[fitresult_c,gof_c] = Fit_Function_c(degree,CC);

% define functions
function [fitresult_p, gof_p] = Fit_Function_p(distr_count,distr_rel)
[xData,yData] = prepareCurveData(distr_count,distr_rel);
ft = fittype('power1');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
[fitresult_p, gof_p] = fit(xData,yData,ft,opts); % goodness of fit

subplot(2,2,4)
plot(fitresult_p,xData,yData);
legend('Data','Power law fit');
xlabel('Degree');
ylabel('Probability');
end
function [fitresult_c, gof_c] = Fit_Function_c(degree,CC)
[xData,yData] = prepareCurveData(degree,CC);
ft = fittype('a*x^(-1)','independent','x','dependent','y');
%excludedPoints = excludedata(xData,yData,'Indices',27);
opts = fitoptions('Method','NonlinearLeastSquares');
%opts.Exclude = excludedPoints;
[fitresult_c,gof_c] = fit(xData,yData,ft,opts);
figure
plot(fitresult_c,xData,yData); %,excludedPoints
legend('Data','1/k fit');
xlabel('degree');
ylabel('clustering coefficient');
end