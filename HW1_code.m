% Clean workspace
clear all; close all; clc
load subdata.mat % Imports the data as the 262144x49 (space by time) matrix called subdata 5
L = 10; % spatial domain n = 64; % Fourier modes x2 = linspace(-L,L,n+1); x = x2(1:n);
y = x;
z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k); [X,Y,Z]=meshgrid(x,y,z); [Kx,Ky,Kz]=meshgrid(ks,ks,ks);
for j=1:49
Un(:,:,:)=reshape(subdata(:,j),n,n,n);
M = max(abs(Un),[],'all');
close all, isosurface(X,Y,Z,abs(Un)/M,0.7) axis([-20 20 -20 20 -20 20]), grid on, drawnow pause(1)
end
ave = zeros(n,n,n); for i = 1:49
ave = ave + fftshift(fftn(reshape(subdata(:,i),n,n,n)));
end
ave = abs(fftshift(ave))/49;
maxAve = max(abs(ave),[],'all');
[x,y,z] = ind2sub(size(ave), find(ave == maxAve)); x_center = Kx(x,y,z)
Appendix B
y_center = Ky(x,y,z) z_center = Kz(x,y,z)
tau = 0.01;
filter = exp((-tau*(Kx - x_center).^2)+(-tau*(Ky - y_center).^2)+(-tau*(Kz - z_center).^2)); dir = zeros(3,49);
for k = 1:49
Un = reshape(subdata(:,k),n,n,n);
unft = filter.*(fftshift(fftn(Un)));
unf = ifftn(unft);
maxF = max(abs(unf),[],'all');
[x,y,z] = ind2sub(size(unf),find(abs(unf) == maxF)); dir(1:k) = Kx(x,y,z);
dir(2:k) = Ky(x,y,z);
dir(3:k) = Kz(x,y,z);
end
subplot(2,2,1) plot3(dir(1,:),dir(2,:),dir(3,:),'-or','Linewidth',2); xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
title('Path of Submarine (tau = 0.01)')
dir(:,end)
tau = 0.1;
filter = exp((-tau*(Kx - x_center).^2)+(-tau*(Ky - y_center).^2)+(-tau*(Kz - z_center).^2)); dir = zeros(3,49);
for k = 1:49
Un = reshape(subdata(:,k),n,n,n);
unft = filter.*(fftshift(fftn(Un)));
unf = ifftn(unft);
maxF = max(abs(unf),[],'all');
[x,y,z] = ind2sub(size(unf),find(abs(unf) == maxF)); dir(1:k) = Kx(x,y,z);
dir(2:k) = Ky(x,y,z);
dir(3:k) = Kz(x,y,z);
end
subplot(2,2,2) plot3(dir(1,:),dir(2,:),dir(3,:),'-or','Linewidth',2); xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
title('Path of Submarine (tau = 0.1)')
dir(:,end)
tau = 1;
filter = exp((-tau*(Kx - x_center).^2)+(-tau*(Ky - y_center).^2)+(-tau*(Kz - z_center).^2)); dir = zeros(3,49);
for k = 1:49
Un = reshape(subdata(:,k),n,n,n);
unft = filter.*(fftshift(fftn(Un)));
unf = ifftn(unft);
maxF = max(abs(unf),[],'all');
[x,y,z] = ind2sub(size(unf),find(abs(unf) == maxF)); dir(1:k) = Kx(x,y,z);

dir(2:k) = Ky(x,y,z);
dir(3:k) = Kz(x,y,z);
end
subplot(2,2,3) plot3(dir(1,:),dir(2,:),dir(3,:),'-or','Linewidth',2); xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
title('Path of Submarine (tau = 1)')
dir(:,end)
tau = 10;
filter = exp((-tau*(Kx - x_center).^2)+(-tau*(Ky - y_center).^2)+(-tau*(Kz - z_center).^2)); dir = zeros(3,49);
for k = 1:49
Un = reshape(subdata(:,k),n,n,n);
unft = filter.*(fftshift(fftn(Un)));
unf = ifftn(unft);
maxF = max(abs(unf),[],'all');
[x,y,z] = ind2sub(size(unf),find(abs(unf) == maxF)); dir(1:k) = Kx(x,y,z);
dir(2:k) = Ky(x,y,z);
dir(3:k) = Kz(x,y,z);
end
subplot(2,2,4) plot3(dir(1,:),dir(2,:),dir(3,:),'-or','Linewidth',2); xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
title('Path of Submarine (tau = 10)')
dir(:,end)