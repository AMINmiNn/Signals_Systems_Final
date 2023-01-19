clc,clear;
A = readmatrix("4_1.csv");
v = A(:,2);
Fs = 6400;  %最高分析至3200Hz,但由于抽样频率较低，可能产生混叠误差，即高频分量叠加到50Hz上。
t = (0:1/Fs:1/Fs*(length(v)-1))';
N0 = 128;  %一个周期约有128个点
r0 = 0.5;
L = 512;  %四个周波窗长

%中间变量
phase = zeros(2,1);
frequency = zeros(2,1);
amptitude = zeros(2,1);
index = zeros(2,1);

%输出变量
FundaAmp = zeros(2,1);
FundaPhase = zeros(2,1);
FundaFrequency = zeros(2,1);


% % 加窗操作
j = 1; %对第i个点求第j次幅值和相位，选取泄露最小的作为i点的结果
for i = 1:(length(v) - N0*10)   %对第i个点求幅值和相位
    for N = N0*9:N0*10   %对每个点加窗，窗长变化以获取最优解
        Xv = fft(v(i:N+i-1).*blackmanharris(N));
        [y2,index2] = max(abs(Xv));  %最大幅值及其位置
        y1 = abs(Xv(index2 + 1));
        k = index2 - 1; %最大幅值对应的k2
        alpha = y2/y1;
        myfun = @(r) deviation(r,alpha);
        r = fzero(myfun,r0);  %偏移量r
        phase(j) = angle(Xv(index2)) - r*pi;
        frequency(j) = (index2 - 1 + r)*Fs/(N);
        amptitude(j)= 2*y2*pi*r*(1-r^2)*(4-r^2)*(9-r^2)/(sin(r*pi)*(12.915-1.22511*r^2 ...
        +0.02913*r^4-0.00006*r^6))/(N);
        j = j + 1;
    end
        [FundaAmp(i),index(i)] = max(amptitude);
        FundaPhase(i) = phase(index(i));
        FundaFrequency(i) = frequency(index(i));
        j = 1;
 end

 for i = (length(v) - N0*10 +1):(length(v) - N0*4)   %对第i个点求幅值和相位,最后两个周波不需要给出
    for N = N0*3:N0*4   %对每个点加窗，窗长变化以获取最优解
        Xv = fft(v(i:N+i-1).*blackmanharris(N));
        [y2,index2] = max(abs(Xv));  %最大幅值及其位置
        y1 = abs(Xv(index2 + 1));
        k = index2 - 1; %最大幅值对应的k2
        alpha = y2/y1;
        myfun = @(r) deviation(r,alpha);
        r = fzero(myfun,r0);  %偏移量r
        phase(j) = angle(Xv(index2)) - r*pi;
        frequency(j) = (index2 - 1 + r)*Fs/(N);
        amptitude(j)= 2*y2*pi*r*(1-r^2)*(4-r^2)*(9-r^2)/(sin(r*pi)*(12.915-1.22511*r^2 ...
        +0.02913*r^4-0.00006*r^6))/(N);
        j = j + 1;
    end
        [FundaAmp(i),index(i)] = max(amptitude);
        FundaPhase(i) = phase(index(i));
        FundaFrequency(i) = frequency(index(i));
        j = 1;
 end

 for i = (length(v) - N0*4 + 1):(length(v) - N0*2)   %对第i个点求幅值和相位,最后两个周波不需要给出
        N = N0*2;   %对每个点加窗，窗长变化以获取最优解
        Xv = fft(v(i:N+i-1).*blackmanharris(N));
        [y2,index2] = max(abs(Xv));  %最大幅值及其位置
        y1 = abs(Xv(index2 + 1));
        k = index2 - 1; %最大幅值对应的k2
        alpha = y2/y1;
        myfun = @(r) deviation(r,alpha);
        r = fzero(myfun,r0);  %偏移量r
        FundaPhase(i) = angle(Xv(index2)) - r*pi;
        FundaFrequency(i) = (index2 - 1 + r)*Fs/(N);
        FundaAmp(i)= 2*y2*pi*r*(1-r^2)*(4-r^2)*(9-r^2)/(sin(r*pi)*(12.915-1.22511*r^2 ...
        +0.02913*r^4-0.00006*r^6))/(N);
 end

amp = FundaAmp.*cos(FundaPhase);
plot(t(1:length(v) - N0*2),amp);
hold on;
plot(t(1:length(v) - N0*2),v(1:length(v) - N0*2));
legend("基波值","实际值");
title("选做三基波幅值、相位");
hold off;

%转化为相对相位
FundaPhase1 = mod(FundaPhase - 2*pi*50*t(1:length(v) - N0*2),2*pi);
FundaPhase1 = FundaPhase1.*(0<=FundaPhase1 & FundaPhase1 <= pi) + (FundaPhase1 - 2*pi).*(pi<FundaPhase1 & FundaPhase1<2*2*pi);   % x in (-pi,pi]
plot(t(1:length(v) - 2*N0),FundaPhase1);
%转化为有效值幅值
FundaAmp1 = FundaAmp/sqrt(2);
A1 = [t(1:length(v) - 2*N0),FundaAmp1,FundaPhase1];

function f = deviation(r,alpha)
f = alpha*(r+3)*(2*r^6-12*r^5-941*r^4+3844*r^3+35041*r^2-77802*r-390632)+...
    (2*r^6-971*r^4+40837*r^2-430500)*(r-4);
end
