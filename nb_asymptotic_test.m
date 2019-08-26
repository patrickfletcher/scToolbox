function p=nb_asymptotic_test(x_a,x_b,size_factor_a,size_factor_b,mu,phi)

alph = size_factor_a * mu / (1 + phi*mu);
beta = (size_factor_b/size_factor_a) * alph;

total = x_a + x_b;

med = median([alph, beta]);
left = ((x_a + 0.5) / total) < med;
right = ~left;

p = ones(length(x_a),1);
p(left)  = 2 *betacdf((x_a(left)  + 0.5)/total(left),alph(left), beta(left));
p(right) = 2 *(1-betacdf((x_a(right) - 0.5)/total(right),alph(right), beta(right)));