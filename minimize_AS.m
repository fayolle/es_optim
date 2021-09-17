% An algorithm of random search
% x_{n+1} = xn  - an*Sn(xn),
% where Sn(xn) differs from the Kiefer-Wolfowitz algorithm in the following way:
% S(xn) = (2*cn)^{-1}*[f(xn+cn*un)-f(xn-cn*un)]*un,
% where un is a random vector of unit length (it is uniformly distributed over the unit sphere).
% Here f(.) is a scalar function to be minimised.
%
% See p. 104 in Alber and Shilman paper
%
% In p. 104 there is a term \eta_n what is it? Random noise? (uniform? normal?)
function [xn,err] = minimize_AS(x0, f, max_iter)
xn = x0;

err = [];
fxn = f(xn);
err = [err; fxn];

dim = size(xn(:), 1);

for i=1:max_iter
    an = 0.2;
    cn = 1.0; %1.0/sqrt(i);
    un = rand_on_sphere(dim);
    XN1 = xn(:) + (cn/2).*un;
    XN2 = xn(:) - (cn/2).*un;
    xn1 = reshape(XN1, size(xn,1), size(xn,2));
    xn2 = reshape(XN2, size(xn,1), size(xn,2));
    g = f(xn1) - f(xn2);
    g = g(:).*un;
    g = reshape(g, size(xn,1), size(xn,2));
    xn = xn - (an/cn)*g;
    
    fxn = f(xn);
    
    errn = fxn;
    err = [err; errn];
end
end


function x = rand_on_sphere(n)
% x uniformly distributed random point on the (n-1)-d unit sphere
r = normrnd(0.0, 1.0, n, 1);
nr = norm(r, 2);
x = r ./ (nr+eps);
end

