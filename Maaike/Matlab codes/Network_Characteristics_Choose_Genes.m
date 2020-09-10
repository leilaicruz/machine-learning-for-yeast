%clear all
close all
clc

% this code extracts which genes are present according to the Repertoire,
% matches this with the interactions and fits the network to the power law
% (builds on Network_Characteristics_Power_Fit)

% import file
file_network = 'Overview interactions.xlsx';   % choose which network you want to use
original_network = xlsread(file_network);      % orginal_network is imported
network = original_network;             % network will be change

% presence of genes
range_p = 'C28:AR28';                   % choose which genes are present in your network
presence = xlsread('Repertoire_Edited.xlsx',range_p);  
[~,size_p]=size(presence);

% make all interactions with absent genes zero
for ii=1:size_p
    if presence(ii)==0      
        network(ii,:)=0;
        network(:,ii)=0;
    end
    ii=ii+1;
end

% delete absent genes
network(~any(network,2),:) = [];    % delete all-zero rows
network(:,~any(network,1)) = [];    % delete all-zero columns

% make network of zeros and ones
binary_network = network./network;
binary_network(isnan(binary_network))=0;

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
distr = zeros(3,max_degree);            % collect all information in one matrix
distr(1,:) = distr_count;                
distr(2,:) = distr_abs;
distr(3,:) = distr_rel;

% plot of degree distribution
figure
subplot(1,3,1)
histogram(degree,'Binwidth',5)          % particular bandwidth or number?
xlabel('Degree')
ylabel('Counts')
subplot(1,3,2)
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

subplot(1,3,3)
plot(fitresult,xData,yData);
legend('Data','Power law fit');
xlabel('Degree');
ylabel('Probability');
end