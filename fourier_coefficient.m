function a = fourier_coefficient(k, n)
    x = linspace(0,1,n);
    midpoints = zeros(1,n-1);
    for j=1:n-1
        midpoints(j) = (x(j)+x(j+1))/2;
    end
    [xx,kk] = meshgrid(midpoints, k);
    vals = potential(midpoints).*exp(2i*pi*kk.*xx);
    a = sum(vals, 2)/n;
end