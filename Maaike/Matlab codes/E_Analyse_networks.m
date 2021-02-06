% This code analyses properties of the whole network.

% INPUT
load('network');    % retrieved from D_Prepare_networks
load('names');      % retrieved from D_Prepare_networks

% OUTPUT
% network_graph
% N_nodes
% N_edges
% avg_degree
% degree_per_node
% N_edges_max
% density
% avg_pathlength
% pathlength
% diameter
% max_degree
% distr_rel
% CC
% avg_CC
% powerfitdata
    % alpha
    % xmin
    % L
    % pvalue
    % gof

network_graph = cell(size(network));
N_nodes = zeros(size(network));
N_edges = zeros(size(network));
avg_degree = zeros(size(network));
degree_per_node = cell(size(network));
N_edges_max = zeros(size(network));
density = zeros(size(network));
avg_pathlength = zeros(size(network));
pathlength = cell(size(network));
diameter = zeros(size(network));
distr_rel = cell(size(network));
CC = cell(size(network));
avg_CC = zeros(size(network));

alpha = zeros(size(network));
xmin = zeros(size(network));
L = zeros(size(network));
pvalue = zeros(size(network));
gof = zeros(size(network));

for ii=1:size(network,2)
    
    % make the network
    network_graph{ii} = graph(network{ii},names{ii});

    % nodes and edges
    N_nodes(ii) = numnodes(network_graph{ii});
    N_edges(ii) = numedges(network_graph{ii});

    % degree
    avg_degree(ii) = 2*N_edges(ii)/N_nodes(ii);                 % average degree <k>=2L/N          
    degree_per_node{ii} = degree(network_graph{ii});

    % density
    N_edges_max(ii) = N_nodes(ii)*(N_nodes(ii)-1)/2;            % max edges = N(N-1)/2
    density(ii) = N_edges(ii)/N_edges_max(ii);
 
    % pathlength
    pathlength{ii} = distances(network_graph{ii});
    avg_pathlength(ii) = sum(pathlength{ii}(~isinf(pathlength{ii})),'all')/nnz(~isinf(pathlength{ii}));
    if ~ isempty(pathlength{ii})
        diameter(ii) = max(max(pathlength{ii}));
    end
      
    % average clustering coefficient calculation (SLOW!)
    CC{ii} = clusteringcoef(network_graph{ii});
    avg_CC(ii) = sum(CC{ii}(~isnan(CC{ii})))/sum(~isnan(CC{ii}));    
 
    % power law fit to degree distribution (SLOW!)
    if N_nodes(ii)>3
        [alpha(ii), xmin(ii), L(ii)] = plfit(degree_per_node{ii});           
        [pvalue(ii),gof(ii)] = plpva(degree_per_node{ii},xmin(ii));
    end
end

all_variables = [avg_CC;avg_degree;avg_pathlength;density;N_edges;N_nodes;alpha;pvalue];

clear ii
clear distr_abs

% define functions
function [alpha, xmin, L] = plfit(x, varargin)
% PLFIT fits a power-law distributional model to data.
%    Source: http://www.santafe.edu/~aaronc/powerlaws/
% 
%    PLFIT(x) estimates x_min and alpha according to the goodness-of-fit
%    based method described in Clauset, Shalizi, Newman (2007). x is a 
%    vector of observations of some quantity to which we wish to fit the 
%    power-law distribution p(x) ~ x^-alpha for x >= xmin.
%    PLFIT automatically detects whether x is composed of real or integer
%    values, and applies the appropriate method. For discrete data, if
%    min(x) > 1000, PLFIT uses the continuous approximation, which is 
%    a reliable in this regime.
%   
%    The fitting procedure works as follows:
%    1) For each possible choice of x_min, we estimate alpha via the 
%       method of maximum likelihood, and calculate the Kolmogorov-Smirnov
%       goodness-of-fit statistic D.
%    2) We then select as our estimate of x_min, the value that gives the
%       minimum value D over all values of x_min.
%
%    Note that this procedure gives no estimate of the uncertainty of the 
%    fitted parameters, nor of the validity of the fit.
%
%    Example:
%       x = (1-rand(10000,1)).^(-1/(2.5-1));
%       [alpha, xmin, L] = plfit(x);
%
%    The output 'alpha' is the maximum likelihood estimate of the scaling
%    exponent, 'xmin' is the estimate of the lower bound of the power-law
%    behavior, and L is the log-likelihood of the data x>=xmin under the
%    fitted power law.
%    
%    For more information, try 'type plfit'
%
%    See also PLVAR, PLPVA

% Version 1.0    (2007 May)
% Version 1.0.2  (2007 September)
% Version 1.0.3  (2007 September)
% Version 1.0.4  (2008 January)
% Version 1.0.5  (2008 March)
% Version 1.0.6  (2008 July)
% Version 1.0.7  (2008 October)
% Version 1.0.8  (2009 February)
% Version 1.0.9  (2009 October)
% Version 1.0.10 (2010 January)
% Version 1.0.11 (2012 January)
% Copyright (C) 2008-2012 Aaron Clauset (Santa Fe Institute)
% Distributed under GPL 2.0
% http://www.gnu.org/copyleft/gpl.html
% PLFIT comes with ABSOLUTELY NO WARRANTY
% 
% Notes:
% 
% 1. In order to implement the integer-based methods in Matlab, the numeric
%    maximization of the log-likelihood function was used. This requires
%    that we specify the range of scaling parameters considered. We set
%    this range to be [1.50 : 0.01 : 3.50] by default. This vector can be
%    set by the user like so,
%    
%       a = plfit(x,'range',[1.001:0.001:5.001]);
%    
% 2. PLFIT can be told to limit the range of values considered as estimates
%    for xmin in three ways. First, it can be instructed to sample these
%    possible values like so,
%    
%       a = plfit(x,'sample',100);
%    
%    which uses 100 uniformly distributed values on the sorted list of
%    unique values in the data set. Second, it can simply omit all
%    candidates above a hard limit, like so
%    
%       a = plfit(x,'limit',3.4);
%    
%    Finally, it can be forced to use a fixed value, like so
%    
%       a = plfit(x,'xmin',3.4);
%    
%    In the case of discrete data, it rounds the limit to the nearest
%    integer.
% 
% 3. When the input sample size is small (e.g., < 100), the continuous 
%    estimator is slightly biased (toward larger values of alpha). To
%    explicitly use an experimental finite-size correction, call PLFIT like
%    so
%    
%       a = plfit(x,'finite');
%    
%    which does a small-size correction to alpha.
%
% 4. For continuous data, PLFIT can return erroneously large estimates of 
%    alpha when xmin is so large that the number of obs x >= xmin is very 
%    small. To prevent this, we can truncate the search over xmin values 
%    before the finite-size bias becomes significant by calling PLFIT as
%    
%       a = plfit(x,'nosmall');
%    
%    which skips values xmin with finite size bias > 0.1.

vec     = [];
sample  = [];
xminx   = [];
limit   = [];
finite  = false;
nosmall = false;
nowarn  = false;

% parse command-line parameters; trap for bad input
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i},
        case 'range',        vec     = varargin{i+1}; i = i + 1;
        case 'sample',       sample  = varargin{i+1}; i = i + 1;
        case 'limit',        limit   = varargin{i+1}; i = i + 1;
        case 'xmin',         xminx   = varargin{i+1}; i = i + 1;
        case 'finite',       finite  = true;
        case 'nowarn',       nowarn  = true;
        case 'nosmall',      nosmall = true;
        otherwise, argok=0; 
    end
  end
  if ~argok, 
    disp(['(PLFIT) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end
if ~isempty(vec) && (~isvector(vec) || min(vec)<=1),
	fprintf('(PLFIT) Error: ''range'' argument must contain a vector; using default.\n');
    vec = [];
end;
if ~isempty(sample) && (~isscalar(sample) || sample<2),
	fprintf('(PLFIT) Error: ''sample'' argument must be a positive integer > 1; using default.\n');
    sample = [];
end;
if ~isempty(limit) && (~isscalar(limit) || limit<min(x)),
	fprintf('(PLFIT) Error: ''limit'' argument must be a positive value >= 1; using default.\n');
    limit = [];
end;
if ~isempty(xminx) && (~isscalar(xminx) || xminx>=max(x)),
	fprintf('(PLFIT) Error: ''xmin'' argument must be a positive value < max(x); using default behavior.\n');
    xminx = [];
end;

% reshape input vector
x = reshape(x,numel(x),1);

% select method (discrete or continuous) for fitting
if     isempty(setdiff(x,floor(x))), f_dattype = 'INTS';
elseif isreal(x),    f_dattype = 'REAL';
else                 f_dattype = 'UNKN';
end;
if strcmp(f_dattype,'INTS') && min(x) > 1000 && length(x)>100,
    f_dattype = 'REAL';
end;

% estimate xmin and alpha, accordingly
switch f_dattype,
    
    case 'REAL',
        xmins = unique(x);
        xmins = xmins(1:end-1);
        if ~isempty(xminx),
            xmins = xmins(find(xmins>=xminx,1,'first'));
        end;
        if ~isempty(limit),
            xmins(xmins>limit) = [];
        end;
        if ~isempty(sample),
            xmins = xmins(unique(round(linspace(1,length(xmins),sample))));
        end;
        dat   = zeros(size(xmins));
        z     = sort(x);
        for xm=1:length(xmins)
            xmin = xmins(xm);
            z    = z(z>=xmin); 
            n    = length(z);
            % estimate alpha using direct MLE
            a    = n ./ sum( log(z./xmin) );
            if nosmall,
                if (a-1)/sqrt(n) > 0.1
                    dat(xm:end) = [];
                    xm = length(xmins)+1;
                    break;
                end;
            end;
            % compute KS statistic
            cx   = (0:n-1)'./n;
            cf   = 1-(xmin./z).^a;
            dat(xm) = max( abs(cf-cx) );
        end;
        D     = min(dat);
        xmin  = xmins(find(dat<=D,1,'first'));
        z     = x(x>=xmin);
        n     = length(z); 
        alpha = 1 + n ./ sum( log(z./xmin) );
        if finite, alpha = alpha*(n-1)/n+1/n; end; % finite-size correction
        if n < 50 && ~finite && ~nowarn,
            fprintf('(PLFIT) Warning: finite-size bias may be present.\n');
        end;
        L = n*log((alpha-1)/xmin) - alpha.*sum(log(z./xmin));

    case 'INTS',
        
        if isempty(vec),
            vec  = (1.50:0.01:3.50);    % covers range of most practical 
        end;                            % scaling parameters
        zvec = zeta(vec);

        xmins = unique(x);
        xmins = xmins(1:end-1);
        if ~isempty(xminx),
            xmins = xmins(find(xmins>=xminx,1,'first'));
        end;
        if ~isempty(limit),
            limit = round(limit);
            xmins(xmins>limit) = [];
        end;
        if ~isempty(sample),
            xmins = xmins(unique(round(linspace(1,length(xmins),sample))));
        end;
        if isempty(xmins)
            fprintf('(PLFIT) Error: x must contain at least two unique values.\n');
            alpha = NaN; xmin = x(1); D = NaN;
            return;
        end;
        xmax   = max(x);
        dat    = zeros(length(xmins),2);
        z      = x;
        fcatch = 0;

        for xm=1:length(xmins)
            xmin = xmins(xm);
            z    = z(z>=xmin);
            n    = length(z);
            % estimate alpha via direct maximization of likelihood function
            if fcatch==0
                try
                    % vectorized version of numerical calculation
                    zdiff = sum( repmat((1:xmin-1)',1,length(vec)).^-repmat(vec,xmin-1,1) ,1);
                    L = -vec.*sum(log(z)) - n.*log(zvec - zdiff);
                catch
                    % catch: force loop to default to iterative version for
                    % remainder of the search
                    fcatch = 1;
                end;
            end;
            if fcatch==1
                % force iterative calculation (more memory efficient, but 
                % can be slower)
                L       = -Inf*ones(size(vec));
                slogz   = sum(log(z));
                xminvec = (1:xmin-1);
                for k=1:length(vec)
                    L(k) = -vec(k)*slogz - n*log(zvec(k) - sum(xminvec.^-vec(k)));
                end
            end;
            [Y,I] = max(L);
            % compute KS statistic
            fit = cumsum((((xmin:xmax).^-vec(I)))./ (zvec(I) - sum((1:xmin-1).^-vec(I))));
            cdi = cumsum(hist(z,xmin:xmax)./n);
            dat(xm,:) = [max(abs( fit - cdi )) vec(I)];
        end
        % select the index for the minimum value of D
        [D,I] = min(dat(:,1));
        xmin  = xmins(I);
        z     = x(x>=xmin);
        n     = length(z);
        alpha = dat(I,2);
        if finite, alpha = alpha*(n-1)/n+1/n; end; % finite-size correction
        if n < 50 && ~finite && ~nowarn,
            fprintf('(PLFIT) Warning: finite-size bias may be present.\n');
        end;
        L     = -alpha*sum(log(z)) - n*log(zvec(find(vec<=alpha,1,'last')) - sum((1:xmin-1).^-alpha));

    otherwise,
        fprintf('(PLFIT) Error: x must contain only reals or only integers.\n');
        alpha = [];
        xmin  = [];
        L     = [];
        return;
end;
end
function [pvalue,gof] = plpva(x, xmin, varargin)
% PLPVA calculates the p-value for the given power-law fit to some data.
%    Source: http://www.santafe.edu/~aaronc/powerlaws/
% 
%    PLPVA(x, xmin) takes data x and given lower cutoff for the power-law
%    behavior xmin and computes the corresponding p-value for the
%    Kolmogorov-Smirnov test, according to the method described in 
%    Clauset, Shalizi, Newman (2007).
%    PLPVA automatically detects whether x is composed of real or integer
%    values, and applies the appropriate method. For discrete data, if
%    min(x) > 1000, PLPVA uses the continuous approximation, which is 
%    a reliable in this regime.
%   
%    The fitting procedure works as follows:
%    1) For each possible choice of x_min, we estimate alpha via the 
%       method of maximum likelihood, and calculate the Kolmogorov-Smirnov
%       goodness-of-fit statistic D.
%    2) We then select as our estimate of x_min, the value that gives the
%       minimum value D over all values of x_min.
%
%    Note that this procedure gives no estimate of the uncertainty of the 
%    fitted parameters, nor of the validity of the fit.
%
%    Example:
%       x = (1-rand(10000,1)).^(-1/(2.5-1));
%       [p, gof] = plpva(x, 1);
%
%    For more information, try 'type plpva'
%
%    See also PLFIT, PLVAR

% Version 1.0   (2007 May)
% Version 1.0.2 (2007 September)
% Version 1.0.3 (2007 September)
% Version 1.0.4 (2008 January)
% Version 1.0.5 (2008 March)
% Version 1.0.6 (2008 April)
% Version 1.0.7 (2009 October)
% Version 1.0.8 (2012 January)
% Copyright (C) 2008-2012 Aaron Clauset (Santa Fe Institute)
% Distributed under GPL 2.0
% http://www.gnu.org/copyleft/gpl.html
% PLPVA comes with ABSOLUTELY NO WARRANTY
% 
% Notes:
% 
% 1. In order to implement the integer-based methods in Matlab, the numeric
%    maximization of the log-likelihood function was used. This requires
%    that we specify the range of scaling parameters considered. We set
%    this range to be [1.50 : 0.01 : 3.50] by default. This vector can be
%    set by the user like so,
%    
%       p = plpva(x, 1,'range',[1.001:0.001:5.001]);
%    
% 2. PLPVA can be told to limit the range of values considered as estimates
%    for xmin in two ways. First, it can be instructed to sample these
%    possible values like so,
%    
%       a = plpva(x,1,'sample',100);
%    
%    which uses 100 uniformly distributed values on the sorted list of
%    unique values in the data set. Second, it can simply omit all
%    candidates above a hard limit, like so
%    
%       a = plpva(x,1,'limit',3.4);
%    
%    Finally, it can be forced to use a fixed value, like so
%    
%       a = plpva(x,1,'xmin',1);
%    
%    In the case of discrete data, it rounds the limit to the nearest
%    integer.
% 
% 3. The default number of semiparametric repetitions of the fitting
% procedure is 1000. This number can be changed like so
%    
%       p = plvar(x, 1,'reps',10000);
% 
% 4. To silence the textual output to the screen, do this
%    
%       p = plpva(x, 1,'reps',10000,'silent');
% 

vec    = [];
sample = [];
limit  = [];
xminx  = [];
Bt     = [];
quiet  = false;
persistent rand_state;

% parse command-line parameters; trap for bad input
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i},
        case 'range',        vec    = varargin{i+1}; i = i + 1;
        case 'sample',       sample = varargin{i+1}; i = i + 1;
        case 'limit',        limit  = varargin{i+1}; i = i + 1;
        case 'xmin',         xminx  = varargin{i+1}; i = i + 1;
        case 'reps',         Bt     = varargin{i+1}; i = i + 1;
        case 'silent',       quiet  = true;
        otherwise, argok=0; 
    end
  end
  if ~argok, 
    disp(['(PLPVA) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end
if ~isempty(vec) && (~isvector(vec) || min(vec)<=1),
	fprintf('(PLPVA) Error: ''range'' argument must contain a vector; using default.\n');
    vec = [];
end;
if ~isempty(sample) && (~isscalar(sample) || sample<2),
	fprintf('(PLPVA) Error: ''sample'' argument must be a positive integer > 1; using default.\n');
    sample = [];
end;
if ~isempty(limit) && (~isscalar(limit) || limit<1),
	fprintf('(PLPVA) Error: ''limit'' argument must be a positive value >= 1; using default.\n');
    limit = [];
end;
if ~isempty(Bt) && (~isscalar(Bt) || Bt<2),
	fprintf('(PLPVA) Error: ''reps'' argument must be a positive value > 1; using default.\n');
    Bt = [];
end;
if ~isempty(xminx) && (~isscalar(xminx) || xminx>=max(x)),
	fprintf('(PLPVA) Error: ''xmin'' argument must be a positive value < max(x); using default behavior.\n');
    xminx = [];
end;

% reshape input vector
x = reshape(x,numel(x),1);

% select method (discrete or continuous) for fitting
if     isempty(setdiff(x,floor(x))), f_dattype = 'INTS';
elseif isreal(x),    f_dattype = 'REAL';
else                 f_dattype = 'UNKN';
end;
if strcmp(f_dattype,'INTS') && min(x) > 1000 && length(x)>100,
    f_dattype = 'REAL';
end;
N = length(x);
x = reshape(x,N,1); % guarantee x is a column vector
if isempty(rand_state)
    rand_state = cputime;
    rand('twister',sum(100*clock));
end;
if isempty(Bt), Bt = 1000; end;
nof = zeros(Bt,1);

if ~quiet,
    fprintf('Power-law Distribution, p-value calculation\n');
    fprintf('   Copyright 2007-2010 Aaron Clauset\n');
    fprintf('   Warning: This can be a slow calculation; please be patient.\n');
    fprintf('   n    = %i\n   xmin = %6.4f\n   reps = %i\n',length(x),xmin,length(nof));
end;
tic;
% estimate xmin and alpha, accordingly
switch f_dattype,
    
    case 'REAL',
        
        % compute D for the empirical distribution
        z     = x(x>=xmin);	nz   = length(z);
        y     = x(x<xmin); 	ny   = length(y);
        alpha = 1 + nz ./ sum( log(z./xmin) );
        cz    = (0:nz-1)'./nz;
        cf    = 1-(xmin./sort(z)).^(alpha-1);
        gof   = max( abs(cz - cf) );
        pz    = nz/N;

        % compute distribution of gofs from semi-parametric bootstrap
        % of entire data set with fit
        for B=1:length(nof)
            % semi-parametric bootstrap of data
            n1 = sum(rand(N,1)>pz);
            q1 = y(ceil(ny.*rand(n1,1)));
            n2 = N-n1;
            q2 = xmin*(1-rand(n2,1)).^(-1/(alpha-1));
            q  = sort([q1; q2]);

            % estimate xmin and alpha via GoF-method
            qmins = unique(q);
            qmins = qmins(1:end-1);
            if ~isempty(xminx),
                qmins = qmins(find(qmins>=xminx,1,'first'));
            end;
            if ~isempty(limit),
                qmins(qmins>limit) = [];
                if isempty(qmins), qmins = min(q); end;
            end;
            if ~isempty(sample),
                qmins = qmins(unique(round(linspace(1,length(qmins),sample))));
            end;
            dat   = zeros(size(qmins));
            for qm=1:length(qmins)
                  qmin = qmins(qm);
                  zq   = q(q>=qmin);
                  nq   = length(zq);
                  a    = nq ./ sum( log(zq./qmin) );
                  cq   = (0:nq-1)'./nq;
                  cf   = 1-(qmin./zq).^a;
                  dat(qm) = max( abs(cq - cf) );
            end;
            if ~quiet,
                fprintf('[%i]\tp = %6.4f\t[%4.2fm]\n',B,sum(nof(1:B)>=gof)./B,toc/60);
            end;
            % store distribution of estimated gof values
            nof(B) = min(dat);
        end;
        pvalue = sum(nof>=gof)./length(nof);

    case 'INTS',

        if isempty(vec),
            vec  = (1.50:0.01:3.50);    % covers range of most practical 
        end;                            % scaling parameters
        zvec = zeta(vec);

        % compute D for the empirical distribution
        z     = x(x>=xmin);	nz   = length(z);	xmax = max(z);
        y     = x(x<xmin); 	ny   = length(y);

        L  = -Inf*ones(size(vec));
        for k=1:length(vec)
            L(k) = -vec(k)*sum(log(z)) - nz*log(zvec(k) - sum((1:xmin-1).^-vec(k)));
        end
        [Y,I] = max(L);
        alpha = vec(I);

        fit = cumsum((((xmin:xmax).^-alpha))./ (zvec(I) - sum((1:xmin-1).^-alpha)));
        cdi = cumsum(hist(z,(xmin:xmax))./nz);
        gof = max(abs( fit - cdi ));
        pz  = nz/N;

        mmax = 20*xmax;
        pdf = [zeros(xmin-1,1); (((xmin:mmax).^-alpha))'./ (zvec(I) - sum((1:xmin-1).^-alpha))];
        cdf = [(1:mmax+1)' [cumsum(pdf); 1]];

        % compute distribution of gofs from semi-parametric bootstrap
        % of entire data set with fit
        for B=1:length(nof)
            % semi-parametric bootstrap of data
            n1 = sum(rand(N,1)>pz);
            q1 = y(ceil(ny.*rand(n1,1)));
            n2 = N-n1;

            % simple discrete zeta generator
            r2 = sort(rand(n2,1));  c = 1;
            q2 = zeros(n2,1);	    k = 1;
            for i=xmin:mmax+1
                while c<=length(r2) && r2(c)<=cdf(i,2), c=c+1; end;
                q2(k:c-1) = i;
                k = c;
                if k>n2, break; end;
            end;
            q = [q1; q2];

            % estimate xmin and alpha via GoF-method
            qmins = unique(q);
            qmins = qmins(1:end-1);
            if ~isempty(xminx),
                qmins = qmins(find(qmins>=xminx,1,'first'));
            end;
            if ~isempty(limit),
                qmins(qmins>limit) = [];
                if isempty(qmins), qmins = min(q); end;
            end;
            if ~isempty(sample),
                qmins = qmins(unique(round(linspace(1,length(qmins),sample))));
            end;
            dat   = zeros(size(qmins));
            qmax  = max(q); zq = q;
            for qm=1:length(qmins)
                qmin = qmins(qm);
                zq   = zq(zq>=qmin);
                nq   = length(zq);
                if nq>1
                    try
                        % vectorized version of numerical calculation
                        zdiff = sum( repmat((1:qmin-1)',1,length(vec)).^-repmat(vec,qmin-1,1) ,1);
                        L = -vec.*sum(log(zq)) - nq.*log(zvec - zdiff);
                    catch
                       % iterative version (more memory efficient, but slower)
                       L       = -Inf*ones(size(vec));
                       slogzq  = sum(log(zq));
                       qminvec = (1:qmin-1);
                       for k=1:length(vec)
                           L(k) = -vec(k)*slogzq - nq*log(zvec(k) - sum(qminvec.^-vec(k)));
                       end;
                    end;
                    [Y,I] = max(L);

                    fit = cumsum((((qmin:qmax).^-vec(I)))./ (zvec(I) - sum((1:qmin-1).^-vec(I))));
                    cdi = cumsum(hist(zq,(qmin:qmax))./nq);
                    dat(qm) = max(abs( fit - cdi ));
                else
                    dat(qm) = -Inf;
                end;

            end
            if ~quiet,
                fprintf('[%i]\tp = %6.4f\t[%4.2fm]\n',B,sum(nof(1:B)>=gof)./B,toc/60);
            end;
            % -- store distribution of estimated gof values
            nof(B) = min(dat);
        end;
        pvalue = sum(nof>=gof)./length(nof);

    otherwise,
        fprintf('(PLPVA) Error: x must contain only reals or only integers.\n');
        pvalue   = [];
        gof = [];
        return;
end;
end
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