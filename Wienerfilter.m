clear; 
clc;
%��ʼ������
L=200;
%�źų���
N=10;
%�˲����Ľ״�
a=0.2;
wn=randn(L,1);
%wn Ϊ���������źŵ������źţ��������һ��L*1�������ɾ����Ԫ��ֵ��%����(0.0,1.0)֮��
sn=zeros(L,1);
%sn Ϊ�źţ�����һ��L*1�������
hn=zeros(N,1);
%hn Ϊϵͳ��λ������Ӧ����һ��N*1�������
rxx=zeros(N,1);
%rxx Ϊ����غ���������һ��N*1 �������
rxd=zeros(N,1);
%rxd Ϊ����غ���������һ��N*1 �������
yn=zeros(L,1);
%yn Ϊ����źţ�����һ��L*1 �������
xt=zeros(L+N,1);
%����һ��(L+N)*1�������
gn=zeros(L,1);
%gn Ϊ yn�� sn��С��������źţ�����һ��L*1 �������%���ݸ�����ʽs(n)=as(n-1)+w(n),���������ź�
for i=2:L
  sn(i,1)=a*sn(i-1,1)+wn(i,1);
end
sn(1,1)=wn(1,1); 
%��ͼ
subplot(3,3,1);
plot(sn,'r');
axis([0,200,-10,10]);
xlabel('ʱ��');
ylabel('����');
title('sn');
%���������źŷ���cd 
cd=(var(wn))/(1-a^2);
%���źż���
x1=awgn(sn,20); 
x2=awgn(sn,10); 
x3=awgn(sn,6);
%��ͼ
subplot(3,3,2);
plot(x1,'g');
axis([0,200,-10,10]);
xlabel('ʱ��');
ylabel('����');
title('x3');
%���������ź��������źŵĻ���غ���,�˴�x1Ϊ�����źţ�snΪ�����ź�
for i=1:N
   for m=i:1:L
       rxd(i,1)=rxd(i,1)+x1(m,1)*sn(m-i+1,1); 
   end
end
%���������źŵ�����غ��� 
for i=1:N
  for m=i:1:L
    rxx(i,1)=rxx(i,1)+x1(m,1)*x1(m-i+1,1); 
  end 
end
%������غ�������������˹���� 
rxx1=toeplitz(rxx);
%���������
irxx=inv(rxx1);
%�����˲���ϵ��h(n)
hn=irxx*rxd;
for i=1:L
  xt(i+N,1)=x1(i,1); 
end
%ʵ���˲� 
for i=1:L
  for m=1:N
      yn(i,1)=yn(i,1)+xt(i+N+1-m,1)*hn(m,1); 
  end
end
%������С��������ź�en 
en=0;
en=cd-(rxd')*hn;
%������С��������ź�gn 
gn=yn-sn;
%�����˲�����ź�ʱ��ͼ 
subplot(3,3,3);
plot(yn);
axis([0,200,-10,10]);
xlabel('ʱ��');
ylabel('����');
title('yn');
%���������ź�������źŶԱ�ͼ 
subplot(3,3,4);
plot(sn,'r');
axis([0,200,-10,10]);
xlabel('ʱ��');
ylabel('����');
title('sn��yn�Ա�'); 
hold on;
plot(yn,'b');
axis([0,200,-10,10]); 
hold off;
%������С��������ź�ͼ figure;
subplot(3,3,5);
plot(gn);
axis([0,200,-2,2]);
xlabel('ʱ��');
ylabel('����');
title('en');
%%�����ƶ�����G
signal_power_in = mean(sn.^2);  
noise_power1 = mean((gn).^2);  
snrin_db = 10 * log10(signal_power_in / noise_power1);  
fprintf('SNR (in dB): %f\n', snrin_db);
signal_power_out = mean(yn.^2);
noise_power2 = mean((gn).^2);
snrout_db = 10 * log10(signal_power_out / noise_power2);  
fprintf('SNR (in dB): %f\n', snrout_db);
G=snrout_db/snrin_db;
fprintf('G: %f\n', G);