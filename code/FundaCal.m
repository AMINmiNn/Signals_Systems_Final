clc,clear;
A = readmatrix("1_1.csv");
t = A(:,1);
v = A(:,2);
Fs = 10e3;
r0 = 0.5;  %偏移量的迭代初值

t0 = 0.4;  %幅值发生阶跃的时刻
FundaFrequence = 49;  %基波频率

%中间变量
amptitude = zeros(2,1);
phase = zeros(2,1);

%输出变量
FundaAmp = zeros(2,1);
FundaPhase = zeros(2,1);

%确定基波初相位
N0 = round(Fs / FundaFrequence);   %一周期的采样点数
M = floor(t0*Fs/N0);
i = 1;

for N = N0*M:t0*Fs
    Xv = fft(v(1:N).*blackmanharris(N));
    [y2,index2] = max(abs(Xv));  %最大幅值及其位置
    y1 = abs(Xv(index2 + 1));
    k = index2 - 1; %最大幅值对应的k
    alpha = y2/y1;
    myfun = @(r) deviation(r,alpha);
    r = fzero(myfun,r0);  %偏移量r
    amptitude(i)= 2*y2*pi*r*(1-r^2)*(4-r^2)*(9-r^2)/(sin(r*pi)*(12.915-1.22511*r^2 ...
        +0.02913*r^4-0.00006*r^6))/(N)/sqrt(2);
    phase(i) = angle(Xv(index2)) - r*pi;
    i = i+1;
end
    [FundaAmp0,index] = max(amptitude);
    FundaPhase0 = phase(index);

L = N0*4;  %设定窗的宽度
%加窗操作
for i = L/2+1:t(end)*Fs-L/2
     v1 = v(i-L/2:i+L/2-1);  %待加窗的信号
     v2 = v1.*blackmanharris(L);   %加窗处理之后
     Xv2 = fft(v2);
    [y2,index2] = max(abs(Xv2));  %最大幅值及其位置
    y1 = abs(Xv2(index2 + 1));
    k = index2 - 1; %最大幅值对应的k
    alpha = y2/y1;
    myfun = @(r) deviation(r,alpha);
    r = fzero(myfun,r0);  %偏移量r    
    FundaAmp(i) = 2*y2*pi*r*(1-r^2)*(4-r^2)*(9-r^2)/(sin(r*pi)*(12.915-1.22511*r^2 ...
        +0.02913*r^4-0.00006*r^6))/L;
end

%幅值修正2————窗既包含阶跃前，又包含阶跃后的部分，阶跃以前的幅值取Amp2(t0*Fs-L/2)，阶跃以后的幅值取Amp2(t0*Fs+L/2)，
FundaAmp(t0*Fs-L/2+1:t0*Fs-1) = FundaAmp(t0*Fs-L/2);
FundaAmp(t0*Fs:t0*Fs+L/2) = FundaAmp(t0*Fs+L/2);
FundaAmp1 = FundaAmp/sqrt(2);


%相位修正2———阶跃附近的相位修正以及取相对相位
FundaPhase = FundaPhase0 + 2*pi*(FundaFrequence-50)*t(1:(end-L/2-1));
FundaPhase1 = mod(FundaPhase,2*pi);
FundaPhase1 = FundaPhase1.*(0<=FundaPhase1 & FundaPhase1 <= pi) + (FundaPhase1 - 2*pi).*(pi<FundaPhase1 & FundaPhase1<2*2*pi);  

subplot(3,1,1)
amp = FundaAmp.*cos(FundaPhase + 2*pi*50*t(1:(end-L/2-1)));
plot(t(L/2+1:length(amp)),amp(L/2+1:end));
hold on;
plot(t(L/2+1:length(amp)),v(L/2+1:length(amp)));
xlabel("时间/s");
ylabel("幅值/V");
legend("基波值","实际值");
title("必做基波幅值、相位");
hold off;

subplot(3,1,2)
plot(t(L/2+1:end-L/2-1),FundaAmp1(L/2+1:end));
title("必做幅值曲线");
xlabel("时间/s");
ylabel("有效值/V");

subplot(3,1,3)
plot(t(L/2+1:end-L/2-1),FundaPhase1(L/2+1:end));
title("必做相位曲线");
xlabel("时间/s");
ylabel("相对相位");

A1 = [t(L/2+1:end-L/2-1),FundaAmp1(L/2+1:end),FundaPhase1(L/2+1:end)];



function f = deviation(r,alpha)
f = alpha*(r+3)*(2*r^6-12*r^5-941*r^4+3844*r^3+35041*r^2-77802*r-390632)+...
    (2*r^6-971*r^4+40837*r^2-430500)*(r-4);
end

