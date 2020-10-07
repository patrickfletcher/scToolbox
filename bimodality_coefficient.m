function BC = bimodality_coefficient(X)

N = size(X,2);

denom_factor = 3*(N-1)^2/((N-2)*(N-3));

skew = skewness(X, 0, 2);
kurt = kurtosis(X, 0, 2) - 3;

BC = (skew.^2 + 1)./(kurt + denom_factor);