function params = outlier_paramset(method)

if ~exist('method','var')
    method='fixed';
end

params.method=method;

switch method
    case 'fixed'
        params.hithr=Inf; %dummy values
        params.lowthr=-Inf;   

    case {'meanstd','geomeanstd'}
        params.ndev_lo=1; 
        params.ndev_hi=1; 
        
    case 'medmad'
        params.ndev_lo=1; 
        params.ndev_hi=1; 
        params.madflag=0;

    case 'prctile'
        params.lowptile=1;
        params.hiptile=1;  

    otherwise
        error('unknown outlier detection method')
end

params.tailoption='both';

params.min_cells=1;
params.logval=false;

params.manual_low.names=[]; 
params.manual_low.vals=[];
params.manual_high.names=[]; 
params.manual_high.vals=[];

params.lowthr_clipval=-Inf;
params.hithr_clipval=Inf;