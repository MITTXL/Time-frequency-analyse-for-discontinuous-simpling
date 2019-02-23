% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 最小二乘替代
% 时间：20180106
% 附属函数脚本：无
% change log：对三种求解方法的分析对比
% (Caution) 求逆运算不可用
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
fs = 512;
% 均匀时间
t = 1/fs:1/fs:2;
N = length(t);
sig = sin(2*pi*(20*t'.^2 + 10*t'));
%% 制造缺失区间
jump_op = 300; % [可调] 
jimp_ed = 320; % [可调] 
t_ls = t([1:jump_op,jimp_ed:end]);
sig_ls = sig([1:jump_op,jimp_ed:end]);
sig_zero = sig;
sig_zero(jump_op+1:jimp_ed-1) = 0;
N_ls = length(t_ls);
continue_rate = 4/4; % [可调] 
%% 无约束最小二乘数值解
freq_range = [fs/N:fs/N:fs/2];
freq_vector = freq_range(1:floor(length(freq_range)*continue_rate));
phi = ones(N_ls, length(freq_vector)+1);
for n = 1:length(freq_vector)
    phi(:,n) = exp(1i*2*pi*freq_vector(n)*t_ls);
end
% 求解LS方法
% 目前只有方法3可用，其他结果近似正太分布，初步判断为复数矩阵求逆问题
% theta = inv(phi'*phi)*phi'*sig_ls;
% theta = phi\sig_ls;
theta = phi'*sig_ls;
%% 重构信号
resample_rete = 100; % [可调] 
t_full = min(t):1/fs/resample_rete:max(t);
freq_range_re = [fs/N/resample_rete:fs/N/resample_rete:fs/2/resample_rete];
phi = ones(length(t_full),length(freq_vector)+1);
for i = 1:length(freq_vector)
    phi(:,i) = exp(1i*2*pi*freq_vector(i)*t_full);
end
sig_re = phi*theta;
real_sig_re = real(sig_re)/max(real(sig_re));
figure,plot(t_full, real_sig_re)
hold on
plot(t,sig,'g.')
plot(t_ls,sig_ls,'r.')
legend('重构信号', '原信号')
hold off
title('EXP LS未加约束')
%% 重构信号与原信号频域对比
ft = abs(fftshift(fft(sig)));
ft_zero = abs(fftshift(fft(sig_zero)));
% 归一化
ft = ft/max(ft);
theta = theta/max(theta);
ft_zero = ft_zero/max(ft_zero);
figure,plot(1:length(theta), abs(theta)),hold on
plot(0:1:(length(ft)/2),ft(length(ft)/2:length(ft))),hold on
plot(0:1:(length(ft_zero)/2),ft_zero(length(ft_zero)/2:length(ft_zero))),hold off
legend('Our Method', 'FFT 未缺失','FFT Zero')


