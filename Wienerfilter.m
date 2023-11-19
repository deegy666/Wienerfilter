clear; 
clc;
%初始化变量
L=200;
%信号长度
N=10;
%滤波器的阶次
a=0.2;
wn=randn(L,1);
%wn 为用于生成信号的噪声信号，随机生成一个L*1矩阵，生成矩阵的元素值在%区间(0.0,1.0)之间
sn=zeros(L,1);
%sn 为信号，生成一个L*1的零矩阵
hn=zeros(N,1);
%hn 为系统单位脉冲响应生成一个N*1的零矩阵
rxx=zeros(N,1);
%rxx 为自相关函数，生成一个N*1 的零矩阵
rxd=zeros(N,1);
%rxd 为互相关函数，生成一个N*1 的零矩阵
yn=zeros(L,1);
%yn 为输出信号，生成一个L*1 的零矩阵
xt=zeros(L+N,1);
%生成一个(L+N)*1的零矩阵
gn=zeros(L,1);
%gn 为 yn与 sn最小距离误差信号，生成一个L*1 的零矩阵%根据给定公式s(n)=as(n-1)+w(n),生成理想信号
for i=2:L
  sn(i,1)=a*sn(i-1,1)+wn(i,1);
end
sn(1,1)=wn(1,1); 
%画图
subplot(3,3,1);
plot(sn,'r');
axis([0,200,-10,10]);
xlabel('时间');
ylabel('幅度');
title('sn');
%生成期望信号方差cd 
cd=(var(wn))/(1-a^2);
%对信号加噪
x1=awgn(sn,20); 
x2=awgn(sn,10); 
x3=awgn(sn,6);
%画图
subplot(3,3,2);
plot(x1,'g');
axis([0,200,-10,10]);
xlabel('时间');
ylabel('幅度');
title('x3');
%生成输入信号与理想信号的互相关函数,此处x1为输入信号，sn为期望信号
for i=1:N
   for m=i:1:L
       rxd(i,1)=rxd(i,1)+x1(m,1)*sn(m-i+1,1); 
   end
end
%生成输入信号的自相关函数 
for i=1:N
  for m=i:1:L
    rxx(i,1)=rxx(i,1)+x1(m,1)*x1(m-i+1,1); 
  end 
end
%将自相关函数生成托普勒斯矩阵 
rxx1=toeplitz(rxx);
%生成逆矩阵
irxx=inv(rxx1);
%生成滤波器系数h(n)
hn=irxx*rxd;
for i=1:L
  xt(i+N,1)=x1(i,1); 
end
%实现滤波 
for i=1:L
  for m=1:N
      yn(i,1)=yn(i,1)+xt(i+N+1-m,1)*hn(m,1); 
  end
end
%计算最小均方误差信号en 
en=0;
en=cd-(rxd')*hn;
%生成最小距离误差信号gn 
gn=yn-sn;
%画出滤波后的信号时域图 
subplot(3,3,3);
plot(yn);
axis([0,200,-10,10]);
xlabel('时间');
ylabel('幅度');
title('yn');
%画出理想信号与输出信号对比图 
subplot(3,3,4);
plot(sn,'r');
axis([0,200,-10,10]);
xlabel('时间');
ylabel('幅度');
title('sn与yn对比'); 
hold on;
plot(yn,'b');
axis([0,200,-10,10]); 
hold off;
%画出最小距离误差信号图 figure;
subplot(3,3,5);
plot(gn);
axis([0,200,-2,2]);
xlabel('时间');
ylabel('幅度');
title('en');
%%计算制度增益G
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