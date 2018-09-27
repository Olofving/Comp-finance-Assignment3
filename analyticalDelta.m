function sol = analyticalDelta(sigma, r, E, T, s)

d1 = ( log(s/E) + (r + 0.5*sigma^2)*T )/(sigma*sqrt(T));

sol = 0.5*(1+erf(d1/sqrt(2)));