function [x, w] = gaussian_quadrature(n)
%gaussian_quadrature Computes the points and weights for Gaussian quadrature
% as defined in Rec. ITU-R P.1411-11 (Section 3.3)

m = floor((n+1) / 2.0);

x = zeros(n, 1);
w = zeros(n, 1);
delta = 1;

for i = 1:m
    x(i) = cos((4*i-1)/(4*n+2) * pi);

    while(1)
        pm1 = 1;
        p = x(i);
        for j = 2:n
            pm2 = pm1;
            pm1 = p;
            p = (2-1.0/j)*x(i)*pm1 - (1-1.0/j)*pm2;
        end
        pp = n*(x(i)*p - pm1)/(x(i).^2-1);
        delta = p/pp;
        x(i) = x(i)-delta;
        if (delta <= eps)
            break
        end
    end

    xlast = x(i) + delta;
    pm1p = (n-1)*(xlast*pm1-pm2)/(xlast.^2 - 1);
    pm1 = pm1 - delta * pm1p;

    w(i) = 2 * (1-x(i).^2)/((n * pm1).^2);
end
for i = m+1 : n
    idx = floor(n/2.0)+m+1-i;
    x(i) = -x(idx);
    w(i) = w(idx);
end

% to make sure x goes from -1 to 1
x = x(end:-1:1);
w = w(end:-1:1);
return 
end


