function [fval, x, err] = minimize_1p1_ES(fitfun, xstart, N, sigma)
% (1+1) ES

% fitfun name or function handle
% N number of iterations
% sigma initial stepsize
% funeval #function evaluations

xn = xstart;
err = [];
fxn = fitfun(xn);
err = [err; fxn];

fval = feval(fitfun, xstart);
x = xstart;
ss=exp(1/3);
ff=exp(1/3 * (-0.27)/(1-0.27));
for i = 1:N
    d = randn(size(xstart));
    f = feval(fitfun, x + sigma*d);
    if f <= fval
        x = x + sigma * d; fval = f; sigma = sigma * ss;
    else
        sigma = sigma * ff;
    end
    
    fx = fitfun(x);
    err = [err; fx];
end
end
