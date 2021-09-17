function [xn,err] = minimize_GLD(x0, f, T, R, r)
% GLD: 
% https://openreview.net/pdf?id=Skep6TVYDB
%

K = log(R/r);
K = floor(K);

xn = x0;

err = [];
fxn = f(xn);
err = [err; fxn];

dim = size(xn(:), 1);

for i=1:T
    vik = zeros(K+1,dim);
    
    for k=1:(K+1)
        rik = 2^(-k+1)*R;
        vik(k,:) = normrnd(0.0, rik, dim, 1);
    end
    
    xn = argmink(f, xn, vik);
    
    fxn = f(xn);
    
    errn = fxn;
    err = [err; errn];
end
end


function x = argmink(f, xt, vik)
% return x = argmin for k of {f(y) | y=xt, y=xt+vik(1), y=xt+vik(2), ...}
x = xt;
fxt = f(xt);
n = size(vik, 1);
for i=1:n
    viki = vik(i,:);
    xi = xt+viki';
    fxi = f(xi);
    if fxi<fxt
       x = xi;
       fxt = fxi;
    end
end
end

