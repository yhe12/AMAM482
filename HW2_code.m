clc;clear all
figure(1)
[y, Fs] = audioread('GNR.m4a');
tr_gnr = length(y)/Fs; % record time in seconds
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Sweet Child O Mine');
p8 = audioplayer(y,Fs); playblocking(p8);
 
L = tr_gnr;
n = length(y);
t2 = linspace(0, L, n + 1);
t = t2(1:n);
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);
a = 1;
tau = 4;
g = exp(- a * (t - tau).^2);
y = y.';
yg = y .* g;
ygt = fft(yg);
ygt_spec = fftshift(abs(ygt));
figure(2)
subplot(3, 1, 1)
plot(t, y, 'k', 'Linewidth', 2)
hold on
plot(t, y, 'm', 'Linewidth', 2)
set(gca, 'Fontsize', 16), xlabel('time (t)'), ylabel('y(t)')
subplot(3, 1, 2)
plot(t, yg, 'k', 'Linewidth', 2)
set(gca, 'Fontsize', 16), xlabel('time (t)'), ylabel('y(t) * g(t -\tau)')
subplot(3, 1, 3)
plot(ks, abs(fftshift(ygt))/max(abs(ygt)), 'r', 'Linewidth', 2); axis([-50 50 0 0.05])
set(gca, 'Fontsize', 16), xlabel('frequency (k)'), ylabel('fft(y(t) * g(t-\tau))')
 
figure(3)
a = [5 1 0.2];
tau = linspace(0, L, 100)
for jj = 1:length(a)
    ygt_spec = [];
    for j =  1:length(tau)
        g =  exp(-a(jj)*(t - tau(j)).^2);
        y = y.';
        yg = y(j) .* g;
        ygt = fft(yg);
        ygt_spec(:, j) = fftshift(abs(ygt));
    end
    subplot(2,2,jj)
    pcolor(tau, ks, log(ygt_spec) + 1)
    shading interp
    set(gca, 'ylim', [-50 50], 'Fontsize', 16)
    colormap(hot)
    xlabel('time (t)'), ylabel('frequency (k)')
    title(['a = ', num2str(a(jj))], 'Fontsize', 16)
end
