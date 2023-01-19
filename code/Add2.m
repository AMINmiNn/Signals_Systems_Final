clc,clear;
A = readmatrix("3_1.csv");
v = A(:,2);
Fs = 10e3;  %最高分析至3200Hz,但由于抽样频率较低，可能产生混叠误差，即高频分量叠加到50Hz上。
t = (0:1/Fs:1/Fs*(length(v)-1))';
r0 = 0.5;  
plot(t,v);

%中间变量
phase = zeros(2,1);
frequency = zeros(2,1);
amptitude = zeros(2,1);

%输出变量
FundaAmp = zeros(2,1);
FundaPhase = zeros(2,1);

%% 求频率并确定基波初相位
N0 =200;   %一周期的采样点数
M = floor(length(t)/N0);
i = 1;

for N = N0*(M-1):length(v)
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
    frequency(i) = (index2 - 1 + r)*Fs/(N);
    i = i+1;
end

    [FundaAmp0,index] = max(amptitude);
    FundaPhase0 = phase(index);
    FundaFrequency = frequency(index);

%% 求基波幅值
N0 = round(Fs/FundaFrequency);
L = N0*4;  %设定窗的宽度
%加窗操作
for i = L/2+1:length(v)-L/2
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


%% 求基波相位
FundaPhase = FundaPhase0 + 2*pi*(FundaFrequency-50)*t(1:(length(v)-L/2));
FundaPhase1 = mod(FundaPhase,2*pi);
FundaPhase1 = FundaPhase1.*(0<=FundaPhase1 & FundaPhase1 <= pi) + (FundaPhase1 - 2*pi).*(pi<FundaPhase1 & FundaPhase1<2*2*pi);  
%% 可视化
subplot(3,1,1)
plot(t(1:(length(v)-L/2)),FundaAmp(1:length(v)-L/2).*cos(FundaPhase1 + 2*pi*50*t(1:(length(v)-L/2))));
hold on;
plot(t,v);
xlabel("时间/s");
ylabel("电压值/V");
legend("基波值","实际值");
title("选做2基波与实际波形比较");
hold off;

subplot(3,1,2)
plot(t(L/2+1:(length(v)-L/2)),FundaAmp(L/2+1:length(v)-L/2)/sqrt(2));
xlabel("时间/s");
ylabel("幅值（有效值）/V");
title("选做2基波幅值-时间变化曲线");

subplot(3,1,3)
plot(t(L/2+1:(length(v)-L/2)),FundaPhase1(L/2+1:end));
xlabel("时间/s");
ylabel("相对相位/rad");
title("选做2基波相对相位-时间变化曲线");

A1 = [t(L/2+1:(length(v)-L/2)),FundaAmp(L/2+1:length(v)-L/2)/sqrt(2),FundaPhase1(L/2+1:end)];

function f = deviation(r,alpha)
f = alpha*(r+3)*(2*r^6-12*r^5-941*r^4+3844*r^3+35041*r^2-77802*r-390632)+...
    (2*r^6-971*r^4+40837*r^2-430500)*(r-4);
end
