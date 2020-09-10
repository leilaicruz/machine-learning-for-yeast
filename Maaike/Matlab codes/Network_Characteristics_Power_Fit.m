clear all
close all
clc

% this code extracts the interactions of the network and fits the network
% to the power law

% import file
title = 'Overview interactions.xlsx';   % choose which network you want to use
network = xlsread(title);               % network is imported

% make matrix of zeros and ones
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
distr = zeros(3,max_degree);
distr(1,:) = distr_count;                
distr(2,:) = distr_abs;
distr(3,:) = distr_rel;

% plot of degree distribution
figure
subplot(1,3,1)
histogram(degree,'Binwidth',5)
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
legend('distr_rel vs. distr_count','Power law fit');
xlabel('Probability');
ylabel('Degree');
end