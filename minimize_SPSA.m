function [xn,err] = minimize_SPSA(x0, f, max_iter)
% not sure what I am doing wrong?
% https://www.jhuapl.edu/SPSA/PDF-SPSA/Matlab-SPSA_Alg.pdf

a = 0.25;
A = 100;
alpha = 1;
c = 1;
gamma = 1.0/6.0;

xn = x0;

err = [];
fxn = f(xn);
err = [err; fxn];

for i=1:max_iter
    %gradphi = f(xn+ys-fxn)-f(xn-ys+fxn);
    
    an = a/(i+A)^alpha;
    cn = c/i^gamma;
    deltan = 2.0.*binornd(1, 0.5, size(x0)) - 1.0;
    %deltan = 2.0*round(rand(size(x0)))-1.0;
    left = f(xn+deltan.*cn);
    right = f(xn-deltan.*cn);
    gradphi = left - right;
    gradphi = gradphi./(2.0*cn.*deltan);
    
    xn = xn - an.*(gradphi);
    
    fxn = f(xn);
    err = [err; fxn];    
end
end


function [theta,err] = minimize_SPSA_2(x0, loss, n)
% https://www.jhuapl.edu/SPSA/PDF-SPSA/Matlab-SPSA_Alg.pdf

a = 0.25;
A = 100;
alpha = 1;
c = 1;
gamma = 1.0/6.0;

p = size(x0,1);

theta = x0;
err = [];
err = [err; loss(theta)];

for k=0:n-1
    ak = a/(k+1+A)^alpha;
    ck = c/(k+1)^gamma;
    delta = 2*round(rand(p,1)) - 1;
    thetaplus = theta + ck*delta;
    thetaminus = theta - ck*delta;
    yplus = loss(thetaplus);
    yminus = loss(thetaminus);
    ghat = (yplus - yminus)./(2*ck*delta);
    theta = theta - ak*ghat;
    
    err = [err; loss(theta)];
end

end

