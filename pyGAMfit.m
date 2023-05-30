function gam = pyGAMfit(X,y,sopts,options)
arguments
    X
    y
    sopts.n_splines=20
    sopts.spline_order=3
    sopts.lam=0.6
    sopts.penalties='auto'
    sopts.basis='ps' %or 'cp' for periodic
    sopts.edge_knots=string(missing) %aka None in python
    options.fit_intercept = true
    options.max_iter=100
    options.tol=1e-4
    options.verbose=false
end
pygam=py.importlib.import_module('pygam.pygam');

% Sline term
%pygam.terms.s(feature, n_splines=20, spline_order=3, lam=0.6,
%penalties='auto', constraints=None, dtype='numerical', basis='ps',
%by=None, edge_knots=None, verbose=False)

sopts.n_splines=int32(sopts.n_splines);
sopts.spline_order=int32(sopts.spline_order);

if isnumeric(sopts.edge_knots)
    sopts.edge_knots=py.numpy.array(sopts.edge_knots,'float32',pyargs('order','C'));
end

s=pygam.SplineTerm(int32(0),...
                   n_splines=sopts.n_splines, ...
                   spline_order=sopts.spline_order, ...
                   lam=sopts.lam, ...
                   basis=sopts.basis, ...
                   edge_knots=sopts.edge_knots);

%classpygam.pygam.LinearGAM(terms='auto', max_iter=100, tol=0.0001,
%scale=None, callbacks=['deviance', 'diffs'], fit_intercept=True,
%verbose=False, **kwargs)

% args = pyargs('terms',s, ...
%     'fit_intercept',options.fit_intercept, ...
%     'max_iter',int32(options.max_iter), ...
%     'tol',options.tol, ...
%     'verbose',options.verbose);

% gam = pygam.LinearGAM(args);

options.max_iter=int32(options.max_iter);
% gam = pygam.LinearGAM(s, options);

gam = pygam.LinearGAM(terms = s, ...
                      fit_intercept=options.fit_intercept, ...
                      max_iter=options.max_iter, ...
                      tol=options.tol, ...
                      verbose=options.verbose);


X = py.numpy.array(X,'float32',pyargs('order','C'));
y = py.numpy.array(y,'float32',pyargs('order','C'));

gam = gam.fit(X, y);