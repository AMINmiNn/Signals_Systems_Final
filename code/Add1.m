%滤波器只能对固定频段的噪声进行去除，无法很好确定幅值、频率的阶跃。
%0.4001s发生频率阶跃
%0.2401s发生幅值阶跃
%阶跃前频率51Hz，阶跃后频率54Hz

clc,clear;
A = readmatrix("2_1.csv");
v = A(:,2);
Fs = 10e3;  %最高分析至3200Hz,但由于抽样频率较低，可能产生混叠误差，即高频分量叠加到50Hz上。
t = (0:1/Fs:1/Fs*(length(v)-1))';

r0 = 0.5;

%% 小波去噪，找到幅值、频率阶跃时刻
%默认前后两周波不出现阶跃
N0 = 200; %一个周波大约200个点
plot(t,v);
hold on;
% v1 = wdenoise(v,"DenoisingMethod","BlockJS");
% plot(t,v1);
% hold on;
% v2 = wdenoise(v,"DenoisingMethod","Bayes");
% plot(t,v2);
% hold on;
v3 = wdenoise(v,"DenoisingMethod","Minimax");
plot(t,v3);
% v4 = wdenoise(v,"DenoisingMethod","UniversalThreshold");
% plot(t,v4);
hold off;
imf = emd(v3);
hht(imf(:,1),Fs);
% [hs, f, t, imfinsf] = hht(imf(:,1),Fs);
% [fmax, index] = max(imfinsf(2*N0:end-2*N0));  %默认前后两周波不出现阶跃
% t1 = t(index+2*N0-1);
% %通过希尔伯特变换获取阶跃时刻
t1 = 0.2401;   %幅值阶跃
t2 = 0.4002;   %频率阶跃


%% 变量
%通过下面的程序获取阶跃前后的基波频率
FundaFrequence1 = 51;  %阶跃前基波频率
FundaFrequence2 = 55;  %阶跃后基波频率
FundaPhase01 = 0.535578910059288;  %阶跃前基波的初始相位
FundaPhase02 = -3.24388678734717; %阶跃后基波的初始相位

% 中间变量
amptitude = zeros(2,1);
phase = zeros(2,1);
frequency = zeros(2,1);

%输出变量
FundaAmp = zeros(2,1);
FundaPhase = zeros(2,1);
%% 求频率阶跃前的频率并确定基波初相位
% N0 = 200;   %一周期的采样点数
% M = floor(t2*Fs/N0);
% i = 1;
% for N = N0*M:t2*Fs
%     Xv = fft(v(1:N).*blackmanharris(N));
%     [y2,index2] = max(abs(Xv));  %最大幅值及其位置
%     y1 = abs(Xv(index2 + 1));
%     k = index2 - 1; %最大幅值对应的k
%     alpha = y2/y1;
%     myfun = @(r) deviation(r,alpha);
%     r = fzero(myfun,r0);  %偏移量r
%     amptitude(i)= 2*y2*pi*r*(1-r^2)*(4-r^2)*(9-r^2)/(sin(r*pi)*(12.915-1.22511*r^2 ...
%         +0.02913*r^4-0.00006*r^6))/(N)/sqrt(2);
%     phase(i) = angle(Xv(index2)) - r*pi;
%     frequency(i) = (index2 - 1 + r)*Fs/(N);
%     i = i+1;
% end
%     [FundaAmp0,index] = max(amptitude);
%     FundaPhase0 = phase(index);
%     FundaFrequency1 = frequency(index);
%% 求频率阶跃后的频率
N0 = 200;
M = floor((length(t)-t2*Fs)/N0);
i = 1;
for N = N0*M:length(t)-t2*Fs
    Xv = fft(v(t2*Fs:N+t2*Fs-1).*blackmanharris(N));
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
    FundaFrequency1 = frequency(index);

%% 求频率阶跃前的基波幅值、相位
N0 = round(Fs/FundaFrequence1);   %一周期的采样点数
L = N0*4;  %设定窗的宽度

%%加窗操作
for i = L/2+1:t2*Fs-L/2
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
FundaAmp(t1*Fs-L/2+1:t1*Fs-1) = FundaAmp(t1*Fs-L/2);
FundaAmp(t1*Fs:t1*Fs+L/2) = FundaAmp(t1*Fs+L/2);
FundaAmp(t2*Fs-L/2+1:t2*Fs+1) = FundaAmp(t2*Fs-L/2);
FundaAmp1 = FundaAmp/sqrt(2);

%相位修正2———阶跃附近的相位修正以及取相对相位
FundaPhase = FundaPhase01 + 2*pi*(FundaFrequence1-50)*t(1:t2*Fs+1);
FundaPhase1 = mod(FundaPhase,2*pi);
FundaPhase1 = FundaPhase1.*(0<=FundaPhase1 & FundaPhase1 <= pi) + (FundaPhase1 - 2*pi).*(pi<FundaPhase1 & FundaPhase1<2*2*pi);  

% subplot(3,1,1)
amp = FundaAmp.*cos(FundaPhase(1:t2*Fs+1) + 2*pi*50*t(1:t2*Fs+1));
% plot(t(L/2+1:length(amp)),amp(L/2+1:length(amp)));
% hold on;
% plot(t(L/2+1:length(amp)),v(L/2+1:length(amp)));
% xlabel("时间/s");
% ylabel("幅值/V");
% legend("基波值","实际值");
% title("必做基波幅值、相位");
% hold off;
% 
% subplot(3,1,2)
% plot(t(L/2+1:t2*Fs+1),FundaAmp1(L/2+1:t2*Fs+1));
% title("必做幅值曲线");
% xlabel("时间/s");
% ylabel("有效值/V");
% 
% subplot(3,1,3)
% plot(t(L/2+1:t2*Fs+1),FundaPhase1(L/2+1:t2*Fs+1));
% title("必做相位曲线");
% xlabel("时间/s");
% ylabel("相对相位");
 
%% 求频率阶跃后的基波幅值、相位
N0 = round(Fs/FundaFrequence2);   %一周期的采样点数
L = N0*4;  %设定窗的宽度

%%加窗操作
for i = t2*Fs+L/2+2:length(t)-L/2
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


%幅值修正
FundaAmp(t2*Fs+2:t2*Fs+L/2+1) = FundaAmp(t2*Fs+L/2+2);
FundaAmp1 = FundaAmp/sqrt(2);

%相位修正2———阶跃附近的相位修正以及取相对相位
FundaPhase(t2*Fs+2:length(t)-L/2) = FundaPhase02 + 2*pi*(FundaFrequence2-50)*t(t2*Fs+2:length(t)-L/2);
FundaPhase1 = mod(FundaPhase,2*pi);
FundaPhase1 = FundaPhase1.*(0<=FundaPhase1 & FundaPhase1 <= pi) + (FundaPhase1 - 2*pi).*(pi<FundaPhase1 & FundaPhase1<2*2*pi);  

subplot(2,2,1)
amp1 = FundaAmp(t2*Fs+2:length(t)-L/2).*cos(FundaPhase(t2*Fs+2:length(t)-L/2) + 2*pi*50*t(t2*Fs+2:length(t)-L/2));
amp = [amp;amp1];
plot(t(1:length(t)-L/2),amp);
hold on;
plot(t(393:length(t)-L/2),v(393:length(t)-L/2));
xlabel("时间/s");
ylabel("幅值/V");
legend("基波值","实际值");
title("选做1基波幅值、相位");
hold off;

subplot(2,2,2)
plot(t(393:length(t)-L/2),FundaAmp1(393:length(t)-L/2));
title("选做1幅值曲线，幅值阶跃时刻：0.2401s");
xlabel("时间/s");
ylabel("有效值/V");

subplot(2,2,3)
plot(t(393:length(t)-L/2),FundaPhase1(393:length(t)-L/2));
title("选做1相位曲线,频率阶跃时刻：0.4002s");
xlabel("时间/s");
ylabel("相对相位");

subplot(2,2,4)
plot(t(393:length(t)-L/2),v(393:length(t)-L/2)/sqrt(2));
hold on;
plot(t(393:length(t)-L/2),FundaAmp1(393:length(t)-L/2).*cos(FundaPhase1(393:length(t)-L/2) + 2*pi*50*t(393:length(t)-L/2)));
hold off;
xlabel("时间/s");
ylabel("幅值（有效值）/V");
legend("基波值","实际值");
title("选做1基波幅值、相位提交数据（有效值）");

A1 = [t(393:length(t)-L/2),FundaAmp1(393:length(t)-L/2),FundaPhase1(393:length(t)-L/2)];

function f = deviation(r,alpha)
f = alpha*(r+3)*(2*r^6-12*r^5-941*r^4+3844*r^3+35041*r^2-77802*r-390632)+...
    (2*r^6-971*r^4+40837*r^2-430500)*(r-4);
end

% 
