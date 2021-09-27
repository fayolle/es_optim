clear all;
close all;
clc;

% Same toy problem as in
m = 2000;
n = 1000;

A = normrnd(0,1,[m, n]);
b = normrnd(0,1,[m, 1]);

x0 = zeros(n, 1);

f = @(x) 0.5/m*norm(A*x - b,2)^2;

% Stochastic variant of Kiefer-Wolfowitz 
[~, err_as] = minimize_AS(x0, f, 1000);

% Vanilla ES 
[~, err_es] = minimize_ES(x0, f, 1000);

% GLD
[~, err_gld] = minimize_GLD(x0, f, 1000, 1, 0.00001);

% RANDOM PURSUIT
[~, ~, funeval, err_rp] = minimize_RP(f, x0, 1000, 0.1);

% 1+1 ES
[~, x, err_1p1] = minimize_1p1_ES(f, x0, 1000, 1.0);

% SPSA 
[~, err_spsa] = minimize_SPSA(x0, f, 1000);

plot([err_as, err_es, err_rp, err_1p1, err_spsa, err_gld]);
legend('AS', 'ES', 'RP', '1+1 ES', 'SPSA', 'GLD');
