clear all
clc
f=1:1:100000;%Ƶ������
w=2*pi*f;
RC=1/(2*pi*15.9*10^3);
for n=1:length(f)
A(n)=abs(1/(1+i*RC*w(n)));%��ֵ˥��ϵ��
P(n)=angle(1/(1+i*RC*w(n))) * 180 / pi;%����ϵ��
end
figure;
subplot(2,1,1);
plot(f,A,'r');%��ֵ����
title('��ֵ˥������');
subplot(2,1,2);
plot(f,P,'blue');%��λ����
title('��λ����');
Func=tf(1,[RC,1]);%ϵͳ�Ĵ��ݺ��� ������Ҫ����G��s��=1/��s^2+2s+1�����Ϳ�����matlab������G=tf��[1],[1 2 1]��;
figure;
bode(Func);%ϵͳ�Ĳ���ͼ
title('����ͼ');
F1=5000;
F2=20000;
F3=100000;
F4=13000;
%
Fs=250000;%������
N=2000;%��������
n=0:N-1;%1��N�е�������
t=n/Fs;%ʱ������
f=n*Fs/N;%Ƶ������
%
%SigIn=sin(2*pi*F1*t);%ԭʼ�ź�
SigIn=sin(2*pi*F1*t)+sin(2*pi*F2*t)+sin(2*pi*F3*t);
%SigIn=square(2*pi*F4*t,50);

[SigOut,Tr] = lsim(Func,SigIn,t);%�˲����ź�ͼ�� 
%The matrix Y has LENGTH(T) rows and as many columns as  outputs in SYS
%����Ҷ�任�����˲�ǰ���ź�
sigInFFT = fft(SigIn,N);
sigInAmp = abs(sigInFFT) * 2 / N;
sigInAmp(1) = sigInAmp(1) / 2;
sigInPhase = angle(sigInFFT) * 180 / pi;

figure;
subplot(3,1,1);
plot(t,SigIn);
title('ԭʼ�ź�ʱ��ͼ��');
subplot(3,1,2);
plot(f,sigInAmp);
title('ԭʼ�ź�Ƶ���ֵ����');
subplot(3,1,3);
plot(f,sigInPhase.*(abs(sigInFFT) >=1e-9));
title('ԭʼ�ź�Ƶ����λ����');

%����Ҷ�任�����˲�����ź�
sigOutFFT = fft(SigOut,N);%FFT�任
sigOutAmp = abs(sigOutFFT) * 2 / N;  
%��FFT����ʱ����ֵ��С��FFTѡ��ĵ����йأ�����Ӱ��������,Ϊ������ʵ�����Ӧ����Ҫ���任��������2����N��
sigOutPhase = angle(sigOutFFT) * 180 / pi;
figure;
subplot(3,1,1);
plot(Tr,SigOut);
title('����ź�ʱ��ͼ��')
subplot(3,1,2);
plot(f,sigOutAmp);
title('����źŷ�ֵ����');
subplot(3,1,3);
plot(f,sigOutPhase.*(abs(sigOutFFT) >=1e-3));
title('����ź���λ����');

%����ؼ���
[a, b] = xcorr(SigIn);
[c, d] = xcorr(SigOut);
figure;
subplot(2,1,1);
plot(b/Fs,a);
title('ԭʼ�ź������ͼ��');
subplot(2,1,2);
plot(d/Fs,c);
title('����ź������ͼ��');

%PSD_WELCH�������㹦�����ܶ�
window=hanning(2000);   %ѡ�õĴ��ڳ���Ҫ��NFFT��N���ĳ�����ͬ
noverlap=32;   %�ֶ������ص��Ĳ������������ȣ�
%[Pxx,fre1]=pwelch(SigIn,N,Fs,window,noverlap,0.95,dflag);   %�����׹���,����0.95�����Ŷȸ����������䣬�޷���ֵ�ǻ��Ƴ���������
%[Pyy,fre2]=pwelch(SigOut,N,Fs,window,noverlap,0.95,dflag); 
NFFT=2000;
[Pxx,fre1] = pwelch(SigIn,window,noverlap,NFFT,Fs);
[Pyy,fre2] = pwelch(SigOut,window,noverlap,NFFT,Fs);
figure;
subplot(2,1,1);
plot(fre1,abs(Pxx)); 
subplot(2,1,2);
plot(fre1,phase(Pxx)); 
%plot(fre1,10*log10(Pxx));  %���ƹ�����
%xlabel('Ƶ��/Hz');ylabel('������/dB');
title('����Welch������ԭʼ�ź�PSDͼ'); grid on
figure;
subplot(2,1,1);
plot(fre2,abs(Pxx)); 
subplot(2,1,2);
plot(fre2,angle(Pxx)); 
%plot(fre2,10*log10(Pyy));  %���ƹ�����
%xlabel('Ƶ��/Hz');ylabel('������/dB');
title('����Welch����������ź�PSDͼ'); grid on



