function [T] = compute_hessian_gaussian(I, halfwindow, sigma)

    n = size(I,1);
    T = ones(n,n,2,2);
    shw = halfwindow*2+1;
    kernelGxx = ones(shw,shw);
    kernelGyy = ones(shw,shw);
    kernelGxy = ones(shw,shw);
    sumxx = 0;
    sumyy = 0;
    sumxy = 0;
    for y = -halfwindow:halfwindow
        for x = -halfwindow:halfwindow
            
            kernelGxx(halfwindow+x+1,halfwindow+y+1) = Gxx(x, y, sigma);
            kernelGyy(halfwindow+x+1,halfwindow+y+1) = Gyy(x, y, sigma);
            kernelGxy(halfwindow+x+1,halfwindow+y+1) = Gxy(x, y, sigma);
            if x~=0 || y~=0
                sumxx = sumxx + kernelGxx(halfwindow+x+1,halfwindow+y+1);
                sumyy = sumyy + kernelGyy(halfwindow+x+1,halfwindow+y+1);
                sumxy = sumxy + kernelGxy(halfwindow+x+1,halfwindow+y+1);
            end
        end
    end
    kernelGxx(halfwindow+1,halfwindow+1) = -sumxx;
    kernelGyy(halfwindow+1,halfwindow+1) = -sumyy;
    kernelGxy(halfwindow+1,halfwindow+1) = -sumxy;
    for x = 1:n
        for y = 1:n
            fxx = 0;
            fyy = 0;
            fxy = 0;
            for dy = -halfwindow:halfwindow
                  for dx = -halfwindow:halfwindow
                      xk = x + dx+1;
                      yk = y + dy+1;
                      if xk < 1
                         xk = -xk+1; 
                      end
                      if xk > n
                          xk = 2*n-xk+1;
                      end
                      if yk < 1
                          yk = -yk+1;
                      end
                      if yk > n
                          yk = 2*n-yk+1;
                      end
                      vk = I(xk,yk);
                      fxx = fxx + kernelGxx(halfwindow-dx+1,halfwindow-dy+1) * vk;
                      fyy = fyy + kernelGyy(halfwindow-dx+1,halfwindow-dy+1) * vk;
                      fxy = fxy + kernelGxy(halfwindow-dx+1,halfwindow-dy+1) * vk;
                  end
            end
            T(x,y,1,1)=fxx;
            T(x,y,1,2)=fxy;
            T(x,y,2,1)=fxy;
            T(x,y,2,2)=fyy;
        end
    end
end

function d = Gxx(x, y, sigma2)
    d = ((x*x-sigma2) / (sigma2*sigma2)) * exp(-(x*x+y*y)/(2*sigma2));
end

function d = Gyy(x, y, sigma2)
    d = ((y*y-sigma2) / (sigma2*sigma2)) * exp(-(x*x+y*y)/(2*sigma2));
end

function d = Gxy(x, y, sigma2)
    d = ((x*y) / (sigma2*sigma2)) * exp(-(x*x+y*y)/(2*sigma2));
end
