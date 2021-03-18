clear all;clc

mtcl = VideoReader('monte_carlo.mov')
skdp = VideoReader('ski_drop.mov')
numFrames = 0
data = cell([],1)
col = zeros(540*960, 379);

while hasFrame(mtcl)
    numFrames = numFrames + 1;
    mtcl2 = readFrame(mtcl);
    data(numFrames) = mtcl2(:,:,3);
    col(:,numFrames) =  reshape(data(numFrames),540*960,3);
end
[U, S, V] = svd(data, "econ")
subplot(2, 1, 1)
plot(diag(S), 'ko', 'Linewidth', 2)
ylabel('\sigma_j')
set(gca, 'Fontsize', 16, 'Xlim', [0.9 delays+0.1])
subplot(2,1,2)
plot(t(1:end-delays+1), V(:,1), 'r', 'Linewidth', 2)
hold on
plot(t(1:end-delays+1),V(:,2),'b--', 'Linewidth', 2)
xlabel('t')
ylabel('v_j(t)')
set(gca,'Fontsize',16)

X1 = V(1:end-1,1:2)';
X2 = V(2:end-1,1:2)';

[U2, S2, V2] = svd(X1, 'econ');
Stilde = U2'*X*V2*diag(1./diag(S2));
[eV, D] = eig(D);
omega = log(mu)/dt;
Phi = U2*eV;
y0 = Phi\X1(:,1);
u_modes = zeros(length(y0), length(t));
for iter = 1:length(t)
    u_modes(:,iter) = y0.*exp(omega*t(iter));
end
u_dmd = real(Phi*u_modes);



while hasFrame(skdp)
    numFrames = numFrames + 1;
    skdp2 = readFrame(skdp);
    data(numFrames) = skdp2(:,:,3);
    col(:,numFrames) =  reshape(data(numFrames),540*960,3);
end
[U, S, V] = svd(data, "econ")
subplot(2, 1, 1)
plot(diag(S), 'ko', 'Linewidth', 2)
ylabel('\sigma_j')
set(gca, 'Fontsize', 16, 'Xlim', [0.9 delays+0.1])
subplot(2,1,2)
plot(t(1:end-delays+1), V(:,1), 'r', 'Linewidth', 2)
hold on
plot(t(1:end-delays+1),V(:,2),'b--', 'Linewidth', 2)
xlabel('t')
ylabel('v_j(t)')
set(gca,'Fontsize',16)

X1 = V(1:end-1,1:2)';
X2 = V(2:end-1,1:2)';

[U2, S2, V2] = svd(X1, 'econ');
Stilde = U2'*X*V2*diag(1./diag(S2));
[eV, D] = eig(D);
omega = log(mu)/dt;
Phi = U2*eV;
y0 = Phi\X1(:,1);
u_modes = zeros(length(y0), length(t));
for iter = 1:length(t)
    u_modes(:,iter) = y0.*exp(omega*t(iter));
end
u_dmd = real(Phi*u_modes);