function model = pyGAMfit(X,y,sopts,modelopts,options)
arguments
    X
    y
    sopts.n_splines=20
    sopts.spline_order=3
    sopts.lam=0.6
    sopts.penalties='auto' % ["derivative", "l2"]
    sopts.basis='ps' %or 'cp' for periodic
    sopts.edge_knots=string(missing) %aka None in python
    modelopts.fit_intercept = true
    modelopts.max_iter=200
    modelopts.tol=1e-4
    modelopts.verbose=false
    options.model = "LinearGAM"
    options.gridLam = []
    options.gridKnots = []
    options.objective = 'auto'
end
pygam=py.importlib.import_module('pygam.pygam');

% ? auto ==   penalties=["derivative", "l2"],
% cellrank: link log dist gamma

% weights? CR uses fate probabilities
% stats? for discovery of traj-assoc. genes. 

%cluster the fits (genes are obs.)

% Sline term
%pygam.terms.s(feature, n_splines=20, spline_order=3, lam=0.6,
%penalties='auto', constraints=None, dtype='numerical', basis='ps',
%by=None, edge_knots=None, verbose=False)

sopts.n_splines=int32(sopts.n_splines);
sopts.spline_order=int32(sopts.spline_order);
sopts.penalties=cellstr(sopts.penalties);
if isnumeric(sopts.edge_knots)
    sopts.edge_knots=py.numpy.array(sopts.edge_knots,'float32',pyargs('order','C'));
end

s=pygam.SplineTerm(int32(0),...
                   n_splines=sopts.n_splines, ...
                   spline_order=sopts.spline_order, ...
                   lam=sopts.lam, ...
                   basis=sopts.basis, ...
                   penalties = sopts.penalties, ...
                   edge_knots=sopts.edge_knots);

%classpygam.pygam.LinearGAM(terms='auto', max_iter=100, tol=0.0001,
%scale=None, callbacks=['deviance', 'diffs'], fit_intercept=True,
%verbose=False, **kwargs)

kwargs = pyargs('terms',s, ...
    'fit_intercept',modelopts.fit_intercept, ...
    'max_iter',int32(modelopts.max_iter), ...
    'tol',modelopts.tol, ...
    'verbose',modelopts.verbose);

switch options.model
    case "GAM"
        model = pygam.GAM(kwargs);
    case "GammaGAM"
        model = pygam.GammaGAM(kwargs);
    case "InvGaussGAM"
    case "LinearGAM"
        model = pygam.LinearGAM(kwargs);
    case "LogisticGAM"
        model = pygam.LogisticGAM(kwargs);
    case "PoissonGAM"
        model = pygam.PoissonGAM(kwargs);
    case "ExpectileGAM"
        model = pygam.ExpectileGAM(kwargs);
end

X = py.numpy.array(X,'float32',pyargs('order','C'));
y = py.numpy.array(y,'float32',pyargs('order','C'));

model = model.fit(X, y);

%disgusting:
doGridL = false;
doGridK = false;
if ~isempty(options.gridLam)
    doGridL = true;
    options.gridLam = py.numpy.array(options.gridLam);
end
if ~isempty(options.gridKnots)
    doGridK = true;
    options.gridKnots = int32(options.gridKnots);
end

if doGridL&&~doGridK
    model.gridsearch(X,y,progress=modelopts.verbose, lam=options.gridLam, objective=options.objective);
end
if ~doGridL&&doGridK
    model.gridsearch(X,y,progress=modelopts.verbose,n_splines=options.gridKnots, objective=options.objective);
end
if doGridL&&doGridK
    model.gridsearch(X,y,progress=modelopts.verbose,lam=options.gridLam, n_splines=options.gridKnots, objective=options.objective);
end
