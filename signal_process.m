clear all
clc
f=1:1:100000;%频率序列
w=2*pi*f;
RC=1/(2*pi*15.9*10^3);
for n=1:length(f)
A(n)=abs(1/(1+i*RC*w(n)));%幅值衰减系数
P(n)=angle(1/(1+i*RC*w(n))) * 180 / pi;%相移系数
end
figure;
subplot(2,1,1);
plot(f,A,'r');%幅值曲线
title('幅值衰减特性');
subplot(2,1,2);
plot(f,P,'blue');%相位曲线
title('相位特性');
Func=tf(1,[RC,1]);%系统的传递函数 比如你要输入G（s）=1/（s^2+2s+1），就可以在matlab中输入G=tf（[1],[1 2 1]）;
figure;
bode(Func);%系统的伯德图
title('伯德图');
F1=5000;
F2=20000;
F3=100000;
F4=13000;
%
Fs=250000;%采样率
N=2000;%采样点数
n=0:N-1;%1行N列的行向量
t=n/Fs;%时间序列
f=n*Fs/N;%频率序列
%
%SigIn=sin(2*pi*F1*t);%原始信号
SigIn=sin(2*pi*F1*t)+sin(2*pi*F2*t)+sin(2*pi*F3*t);
%SigIn=square(2*pi*F4*t,50);

[SigOut,Tr] = lsim(Func,SigIn,t);%滤波后信号图像 
%The matrix Y has LENGTH(T) rows and as many columns as  outputs in SYS
%傅里叶变换分析滤波前的信号
sigInFFT = fft(SigIn,N);
sigInAmp = abs(sigInFFT) * 2 / N;
sigInAmp(1) = sigInAmp(1) / 2;
sigInPhase = angle(sigInFFT) * 180 / pi;

figure;
subplot(3,1,1);
plot(t,SigIn);
title('原始信号时域图像');
subplot(3,1,2);
plot(f,sigInAmp);
title('原始信号频域幅值特性');
subplot(3,1,3);
plot(f,sigInPhase.*(abs(sigInFFT) >=1e-9));
title('原始信号频域相位特性');

%傅里叶变换分析滤波后的信号
sigOutFFT = fft(SigOut,N);%FFT变换
sigOutAmp = abs(sigOutFFT) * 2 / N;  
%做FFT分析时，幅值大小与FFT选择的点数有关，但不影响分析结果,为了与真实振幅对应，需要将变换后结果乘以2除以N。
sigOutPhase = angle(sigOutFFT) * 180 / pi;
figure;
subplot(3,1,1);
plot(Tr,SigOut);
title('输出信号时域图像')
subplot(3,1,2);
plot(f,sigOutAmp);
title('输出信号幅值特性');
subplot(3,1,3);
plot(f,sigOutPhase.*(abs(sigOutFFT) >=1e-3));
title('输出信号相位特性');

%自相关计算
[a, b] = xcorr(SigIn);
[c, d] = xcorr(SigOut);
figure;
subplot(2,1,1);
plot(b/Fs,a);
title('原始信号自相关图像');
subplot(2,1,2);
plot(d/Fs,c);
title('输出信号自相关图像');

%PSD_WELCH方法计算功率谱密度
window=hanning(2000);   %选用的窗口长度要与NFFT（N）的长度相同
noverlap=32;   %分段序列重叠的采样点数（长度）
%[Pxx,fre1]=pwelch(SigIn,N,Fs,window,noverlap,0.95,dflag);   %功率谱估计,并以0.95的置信度给出置信区间，无返回值是绘制出置信区间
%[Pyy,fre2]=pwelch(SigOut,N,Fs,window,noverlap,0.95,dflag); 
NFFT=2000;
[Pxx,fre1] = pwelch(SigIn,window,noverlap,NFFT,Fs);
[Pyy,fre2] = pwelch(SigOut,window,noverlap,NFFT,Fs);
figure;
subplot(2,1,1);
plot(fre1,abs(Pxx)); 
subplot(2,1,2);
plot(fre1,phase(Pxx)); 
%plot(fre1,10*log10(Pxx));  %绘制功率谱
%xlabel('频率/Hz');ylabel('功率谱/dB');
title('基于Welch方法的原始信号PSD图'); grid on
figure;
subplot(2,1,1);
plot(fre2,abs(Pxx)); 
subplot(2,1,2);
plot(fre2,angle(Pxx)); 
%plot(fre2,10*log10(Pyy));  %绘制功率谱
%xlabel('频率/Hz');ylabel('功率谱/dB');
title('基于Welch方法的输出信号PSD图'); grid on



