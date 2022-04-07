function BC = bimodality_coefficient(X, dim)

N = size(X,dim);

denom_factor = 3*(N-1)^2/((N-2)*(N-3));

skew = skewness(X, 0, dim);
kurt = kurtosis(X, 0, dim) - 3;

BC = (skew.^2 + 1)./(kurt + denom_factor);