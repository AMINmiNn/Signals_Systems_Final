clf;
clear;
% 分成阶跃前和阶跃后两部分进行分析
dataT = readmatrix('1_1.csv');
dataT = dataT.';
fs = 10e3;  % 采样频率
[peakIndex, part1, part2] = getPeakIndex(dataT);  % 拆分出两部分，只取整周期
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

[HarmResult, FundaResult] = anaHarm(fftWinDataFirst)
[HarmResult, FundaResult] = anaHarm(fftWinDataSecond)



function [HarmResult, FundaResult] = anaHarm(fftWinData)
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
    figure(8);
    stem(xIndex, value(1, :));
    hold on;
    stem(xIndex, value(2, :), 'or');
    hold on;
    stem(xIndex, value(3, :));
    
    % 设定阈值为2，最大幅值小于2的一律认为是误判，删去误判的谐波
    level = 2;
    delX = find(value(2,:)<level);
    value(:, delX) = [];  % 删除
    harmIndex(:, delX) = [];  % 删除 
    
    % 重新画图，看各个谐波分量及其邻居，此时没有误判的，而且剩余的谐波数量很少
    xIndex = linspace(1, length(value), length(value));  % 横坐标
    figure(9);
    stem(xIndex, value(1, :));
    hold on;
    stem(xIndex, value(2, :), 'or');
    hold on;
    stem(xIndex, value(3, :));
    
    % 修正谐波频率及幅值，按照blackmanharris窗的修正方程计算
    [delta, vDelta] = blackmanHarrisFix(value);  % 计算修正量
    fixHarmFreq = (harmIndex(2,:)+delta - 1)./dataN.*fs;  % 这里要减1，因为fft结果角标应该从0开始
    fixValue = (value(1,:)+2.*value(2,:)+value(3,:)).*vDelta./dataN;

    % 修正基波频率及幅值
    [fDelta, fVDelta] = blackmanHarrisFix(RealfundaValue);
    fixFundaFreq = (RealfundaIndex(2, :)+fDelta - 1)./dataN.*fs;  % 这里要减一，因为fft结果角标应该从0开始
    fixFundaValue = (RealfundaValue(1,:)+2.*RealfundaValue(2,:)+RealfundaValue(3,:)).*fVDelta./dataN;

    format longE;
    HarmResult = [fixHarmFreq; fixValue];
    FundaResult = [fixFundaFreq; fixFundaValue];
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
function [peakIndex, part1, part2] = getPeakIndex(dataT)
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
        if dataT(2, peakIndex(i)) > dataT(2, peakIndex(i-1)).*1.1
            cutIndex = i;
            break;
        end
    end
    part1 = peakIndex(1:cutIndex-1);
    part2 = peakIndex(cutIndex:end);
end