function [xn,err] = minimize_ES(x0, f, max_iter)
xn = x0;

err = [];
fxn = f(xn);
err = [err; fxn];

sigma = 0.1;
beta = 1.0;
lr = 0.2;
n = size(x0, 1);
scale = sigma / sqrt(n);

for i=1:max_iter
    epsilon = normrnd(0, scale, [n, 1]);
    
    XN1 = xn(:) + epsilon;
    XN2 = xn(:) - epsilon;
    xn1 = reshape(XN1, size(xn,1), size(xn,2));
    xn2 = reshape(XN2, size(xn,1), size(xn,2));
    g = f(xn1) - f(xn2);
    g = g(:).*epsilon;
    g = (beta/(2*sigma^2)).*g;    
    g = reshape(g, size(xn,1), size(xn,2));
    xn = xn - lr.*g;
    
    fxn = f(xn);
    
    errn = fxn;
    err = [err; errn];
end
end
