%必做任务

function [Phasor_funda,Phasor_harmo]=Phasor_Calculation(data)
% data是原始输入信号,2列矩阵,第一列是时间,第二列是电压瞬时值
%Phasor_funda是基波相量计算结果，是一个3列的矩阵，第一列时间，第二列幅值，第三列相位
%Phasor_harmo是谐波计算结果，是一个2列矩阵，第一列是频率，第二列是幅值
t = data(:,1);
v = data(:,2);

%抽样频率
Fs = 10e3;

%偏移量的迭代初值
r0 = 0.5;  

%一周期的采样点数
N0 = 200;   

%找幅值阶跃时刻
t0 =  amptitude_step_find(v,N0,Fs); 

%谐波分析
[Phasor_harmo, Funda] = anaHarmAllFromData(data,t0,Fs);

%中间变量
FundaFrequency = Funda(1);
amptitude = zeros(2,1);
phase = zeros(2,1);

%输出变量
FundaAmp = zeros(2,1);
FundaPhase = zeros(2,1);

%确定基波初相位
N0 = round(Fs/FundaFrequency);
M = floor(t0*Fs/N0);
i = 1;
j = 1;

%通过循环来确定泄露时的抽样长度，进而找到最精确的基波初相位
for N = N0*M:t0*Fs
    Xv = fft(v(1:N).*blackmanharris(N));
    [y2,index2] = max(abs(Xv));  %最大幅值及其位置
    y1 = abs(Xv(index2 + 1));
    alpha = y2/y1;
    myfun = @(r) deviation(r,alpha);
    r = fzero(myfun,r0);  %偏移量r
    amptitude(j)= 2*y2*pi*r*(1-r^2)*(4-r^2)*(9-r^2)/(sin(r*pi)*(12.915-1.22511*r^2 ...
        +0.02913*r^4-0.00006*r^6))/(N)/sqrt(2);
    phase(j) = angle(Xv(index2)) - r*pi;
    j = j+1;
end
    [~,index0] = max(amptitude);
    FundaPhase0 = phase(index0);

%加窗找每个点的幅值
L = N0*4;  %设定窗的宽度
for  j = L/2+1:t(end)*Fs-L/2
     v1 = v(j-L/2:j+L/2-1);  %待加窗的信号
     v2 = v1.*blackmanharris(L);   %加窗处理之后
     Xv2 = fft(v2);
    [y2,index2] = max(abs(Xv2));  %最大幅值及其位置
    y1 = abs(Xv2(index2 + 1));
    alpha = y2/y1;
    myfun = @(r) deviation(r,alpha);
    r = fzero(myfun,r0);  %偏移量r    
    FundaAmp(j) = 2*y2*pi*r*(1-r^2)*(4-r^2)*(9-r^2)/(sin(r*pi)*(12.915-1.22511*r^2 ...
        +0.02913*r^4-0.00006*r^6))/L;
end

%幅值修正2————窗既包含阶跃前，又包含阶跃后的部分，阶跃以前的幅值取Amp2(t0*Fs-L/2)，阶跃以后的幅值取Amp2(t0*+L/2)，
FundaAmp(t0*Fs-L/2+1:t0*Fs-1) = FundaAmp(t0*Fs-L/2);
FundaAmp(t0*Fs:t0*Fs+L/2) = FundaAmp(t0*Fs+L/2);
FundaAmp1 = FundaAmp/sqrt(2);


%相位修正2———阶跃附近的相位修正以及取相对相位
FundaPhase = FundaPhase0 + 2*pi*(FundaFrequency-50)*t(1:(end-L/2-1));
FundaPhase1 = mod(FundaPhase,2*pi);
FundaPhase1 = FundaPhase1.*(0<=FundaPhase1 & FundaPhase1 <= pi) + (FundaPhase1 - 2*pi).*(pi<FundaPhase1 & FundaPhase1<2*2*pi);  

figure(1)
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

Phasor_funda = [t(L/2+1:end-L/2-1),FundaAmp1(L/2+1:end),FundaPhase1(L/2+1:end)];

    function f = deviation(r,alpha)
    f = alpha*(r+3)*(2*r^6-12*r^5-941*r^4+3844*r^3+35041*r^2-77802*r-390632)+...
        (2*r^6-971*r^4+40837*r^2-430500)*(r-4);
    end

    function t0 = amptitude_step_find(v,N0,Fs)
        imf = emd(v);%分离出固有模态函数
        [~, ~, t1, imfinsf] = hht(imf(:,1),Fs);
        [~, index] = max(imfinsf(2*N0:end-2*N0));  %默认前后两周波不出现阶跃
        t0 = t1(index + 2*N0 - 1); %通过希尔伯特变换获取阶跃时刻
    end

    
    function [Harm, Funda] = anaHarmAllFromData(data,t0,Fs)
    dataT = data.';
    fs = 10e3;  % 采样频率
    [peakIndex, part1, part2] = getPeakIndex(dataT,t0,Fs);  % 拆分出两部分，只取整周期
    dataFirst = dataT(:, part1(1):part1(end));
    dataSecond = dataT(:, part2(1):part2(end));

    % 计算序列长度
    nFirst = length(dataFirst);  
    nSecond = length(dataSecond);
    % 计算blackmanharris窗
    winFirst = blackmanharris(nFirst, 'periodic');
    winSecond = blackmanharris(nSecond, 'periodic');
    % 计算加窗后的时域数据
    winDataFirst = dataFirst(2, :).*winFirst';
    winDataSecond = dataSecond(2, :).*winSecond';
    

    % 对加窗后的结果进行fft
    fftWinDataFirst = fft(winDataFirst);
    fftWinDataSecond = fft(winDataSecond);
    

    % 计算原始的fft
    fftDataFirst = fft(dataFirst(2, :));
    fftDataSecond = fft(dataSecond(2, :));

    % 阶跃前的分析结果
    [HarmResult1, FundaResult1] = anaHarm(fftWinDataFirst, fftDataFirst);
    % 阶跃后的分析结果
    %[HarmResult2, FundaResult2] = anaHarm(fftWinDataSecond, fftDataSecond);
    % 将阶跃前的频率，幅值写入文件
   
    Harm = HarmResult1(2:3, :);
    Harm = [Harm(1, :); Harm(2, :)./sqrt(2)];
    Funda = FundaResult1(2:3, :);
    Funda = [Funda(1, :); Funda(2, :)./sqrt(2)];
    end


    function [HarmResult, FundaResult] = anaHarm(fftWinData, fftData)
        fs = 10e3;  % 采样频率
        dataN = length(fftWinData);  % 获取长度
        harmIndex = getHarmIndex(fftWinData);  % 获取谐波角标（此时的谐波仍然包含基波）
        value = abs(fftWinData(harmIndex));  % 获取谐波幅值
        [fundaValue, fundaIndex] = max(sum(value, 1), [], 2);  % 取出基波角标
        RealfundaIndex = harmIndex(:, fundaIndex);  % 获取基波在fft中的角标
        RealfundaValue = value(:, fundaIndex);  % 获取基波幅值
        value(:, fundaIndex) = [];  % 删除基波幅值
        harmIndex(:, fundaIndex) = [];  % 删除基波角标
        
        % 绘图，把谐波分量及邻居都画出来，会发现有误判的，有的分量幅值很小
        xIndex = linspace(1, length(value), length(value));
        figure(2);
        stem(xIndex, value(1, :));
        hold on;
        stem(xIndex, value(2, :), 'or');
        hold on;
        stem(xIndex, value(3, :));
        title('未消除误判前的谐波')
        
        % 设定阈值为2，最大幅值小于2的一律认为是误判，删去误判的谐波
        level = 2;
        delX = find(value(2,:)<level);
        value(:, delX) = [];  % 删除
        harmIndex(:, delX) = [];  % 删除 
        
        % 重新画图，看各个谐波分量及其邻居，此时没有误判的，而且剩余的谐波数量很少
        xIndex = linspace(1, length(value), length(value));  % 横坐标
        figure(3);
        stem(xIndex, value(1, :));
        hold on;
        stem(xIndex, value(2, :), 'or');
        hold on;
        stem(xIndex, value(3, :));
        title('消除误判后的谐波')
        
        % 修正谐波频率及幅值，按照blackmanharris窗的修正方程计算
        [delta, vDelta] = blackmanHarrisFix(value);  % 计算修正量
        fixHarmFreq = (harmIndex(2,:)+delta - 1)./dataN.*fs;  % 这里要减1，因为fft结果角标应该从0开始
        fixValue = (value(1,:)+2.*value(2,:)+value(3,:)).*vDelta./dataN;
    
        % 修正基波频率及幅值
        [fDelta, fVDelta] = blackmanHarrisFix(RealfundaValue);
        fixFundaFreq = (RealfundaIndex(2, :)+fDelta - 1)./dataN.*fs;  % 这里要减一，因为fft结果角标应该从0开始
        fixFundaValue = (RealfundaValue(1,:)+2.*RealfundaValue(2,:)+RealfundaValue(3,:)).*fVDelta./dataN;
        
        % 获取原始相位
        rawAngle = angle(fftData);
        % 计算修正后的谐波相位
        harmAngle = rawAngle(harmIndex(2, :))  - delta.*pi;
        % 计算修正后的基波相位
        fundaAngle = rawAngle(RealfundaIndex(2, :))  - fDelta.*pi;
    
    
        format longE;
        HarmResult = [harmIndex(2, :); fixHarmFreq; fixValue; harmAngle];  % 原始角标 修正频率 修正幅值 修正相位
        FundaResult = [RealfundaIndex(2, :); fixFundaFreq; fixFundaValue; fundaAngle];
    end
    
    
    % blackmanHarris窗三次插值修正函数
    function [delta, vDelta] = blackmanHarrisFix(y)
        alpha = (y(3,:)-y(1,:))./y(2,:);
        delta = 0.93891885.*alpha - 0.08203810.*alpha.^3 + 0.01541147.*alpha.^5 - 0.00316988.*alpha.^7;
        vDelta = 1.65866360 + 0.44865377.*delta.^2 + 0.06478728.*delta.^4 + 0.00693833.*delta.^6;
    end
    
    
    % 找最大幅值及其左右两个次大幅值的角标
    function harmIndex = getHarmIndex(fftData)
        fftData = abs(fftData);
        len = length(fftData)/2;  % 获取长度，取一半即可，因为是对称的
        harmIndex = zeros(3, 5);  % 第一行为左邻居，第三行为右邻居
        index = 1;  % 存储用的角标
        for i = 1 + 1: len - 1
            if fftData(i) > fftData(i-1) && fftData(i) > fftData(i+1)
                harmIndex(1, index) = i-1;
                harmIndex(2, index) = i;
                harmIndex(3, index) = i+1;
                index = index + 1;  % 角标自增
            end
        end
    end
    
    
    % 找到左右峰值对应的角标
    function [peakIndex, part1, part2] = getPeakIndex(dataT,t0,Fs)
        %peakIndex 返回了所有峰值的角标，part1返回了所有未阶跃的峰值的角标，part2返回了所有阶跃后的峰值的角标
        len = length(dataT(1, :));  % 获取时域长度
        peakIndex = zeros(1, 10);  % 存储峰值角标
        saveIndex = 1;  % 辅助存储的角标
        margin = 10;  % 观察域，峰值要比左右各margin个点大
        for i = 1+margin: len-margin
            if dataT(2, i) > max(dataT(2, i-margin: i-1)) && dataT(2, i) > max(dataT(2, i+1:i+margin))
                peakIndex(saveIndex) = i;
                saveIndex = saveIndex + 1;   % 角标自增1
            end
        end
        
        cutIndex = 0;  % 分隔阶跃前和阶跃后的角标
        for i = 2:length(peakIndex)  % 若前一个峰值的1.1倍小于后一个峰值，则两个峰值显著不同，在此处进行划分
            if peakIndex(i) >= t0*Fs
                cutIndex = i;
                break;
            end
        end
        part1 = peakIndex(1:cutIndex-1);
        part2 = peakIndex(cutIndex:end);
    end

    writematrix([Phasor_funda,[Phasor_harmo';zeros(size(Phasor_funda,1)-size(Phasor_harmo,2),2)]],'1_1泛化验证.csv');
end
