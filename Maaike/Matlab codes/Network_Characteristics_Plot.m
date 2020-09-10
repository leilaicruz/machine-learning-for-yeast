clear all
close all
clc

% this code extracts which genes are present according to the Repertoire,
% matches this with the interactions and fits the network to the power law
% and visualises it
% (builds on Network_Characteristics_Choose_Genes)

% import file of network
file_network = 'Overview interactions.xlsx';       % choose which network you want to use
original_network = xlsread(file_network);          % orginal_network is imported
network = original_network;                        % network will be change
[~,gene_names,~] = xlsread(file_network,'B1:AQ1'); % string of names of genes

% presence of genes
file_presence ='Repertoire_Edited.xlsx';
range_p = 'C16:AR16';                              % choose which genes are present in your network
presence = xlsread(file_presence,range_p);
[~,size_p]=size(presence);

% delete all interactions with absent genes
for ii=size_p:-1:1                                 % backwards, otherwise for-loop stagnates
    if presence(ii)==0      
        network(ii,:)=[];
        network(:,ii)=[];
        gene_names(ii)=[];
    end
end

% make matrix of zeros and ones
binary_network = network./network;
binary_network(isnan(binary_network))=0;

% visualise the network
visualisation = graph(binary_network,gene_names);
figure
subplot(2,2,1)
plot(visualisation)

% network measures
[N_nodes,~] = size(network);    % number of nodes is number of rows or columns
N_links = nnz(network);         % number of links is number of non-zero elements
avg_degree = 2*N_links/N_nodes; % average degree <k>=2L/N
degree = sum(binary_network);   % number of links per node

%degree distribution
max_degree = max(degree);
distr_count = [1:max_degree];           % nodes with no connections are not taken into account
distr_abs = histcounts(degree);         % absolute distribution
distr_rel = histcounts(degree)/N_nodes; % relative distribution
distr = zeros(3,max_degree);            % create an overview table
distr(1,:) = distr_count;                
distr(2,:) = distr_abs;
distr(3,:) = distr_rel;

% plot of degree distribution
subplot(2,2,2)
histogram(degree,'Binwidth',5)
xlabel('Degree')
ylabel('Counts')
subplot(2,2,3)
scatter(distr_count,distr_rel,'x')
xlabel('Degree')
ylabel('Probability')

% power law fit to degree distribution
[fitresult,gof] = Fit_Function(distr_count,distr_rel);
function [fitresult, gof] = Fit_Function(distr_count,distr_rel)
[xData,yData] = prepareCurveData(distr_count,distr_rel);
ft = fittype('power1');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
[fitresult, gof] = fit(xData,yData,ft,opts); % goodness of fit

subplot(2,2,4)
plot(fitresult,xData,yData);
legend('Data','Power law fit');
xlabel('Degree');
ylabel('Probability');
end