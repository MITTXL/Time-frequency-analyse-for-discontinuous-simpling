% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 最小二乘替代FFT
% 时间：20180106
% 附属函数脚本：无
% change log：非均匀时间序列
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear
fs = 800;
% 均匀时间
t_even = 0:1/fs:1-1/fs;
t_un = sort(rand(1,fs));


y_even = sin(2*pi*(30*t_even'.^2+40*t_even'));
y_un = sin(2*pi*(30*t_un'.^2+40*t_un'));

t = t_un;
y = y_un;
N = length(t);


phi = ones(length(t),N-1);
for n = 1:N/2-1
    phi(:,2*n-1) = sin(2*pi*n*fs/N*t);
    phi(:,2*n) = cos(2*pi*n*fs/N*t);
end

% 核查theta公式是否错误
% theta = (phi'*phi)\phi'*y_un;
theta = phi'*y;
y_hat = phi*theta;


ft = abs(fftshift(fft(y)));
L = zeros(N/2-1,1);
L(1) = theta(end);
for n = 1:N/2-1
    L(n+1) = sqrt(theta(2*n-1)^2 + theta(2*n)^2);
end

ft = ft/max(ft);
L = L/max(L);
figure,plot(L),hold on
plot(0:1:(length(ft)/2),ft(length(ft)/2:length(ft))),hold off

err = abs(y - y_hat);
