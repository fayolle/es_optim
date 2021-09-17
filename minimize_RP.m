function [fval, x, funeval, err] = minimize_RP(fitfun, xstart, N, mu)
% Algorithm Random Pursuit

% fitfun name or function handle
% N number of iterations
% mu line search accuracy
% funeval #function evaluations

xn = xstart;
err = [];
fxn = fitfun(xn);
err = [err; fxn];


%line search parameters
opts = optimset('Display', 'off', 'LargeScale', 'off', 'TolX', mu);
funeval = 0;
x = xstart;

for i = 1:N
    d=randn(size(xstart)); d=d/norm(d);
    [sigma, fval, ~, infos] = ...
        fminunc(@(sigma) feval(fitfun, x + sigma*d), 0, opts);
    funeval = funeval + infos.funcCount;
    x = x + sigma*d;
    
    err = [err; fval];
    
end
end
