clf;

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

% 构造fft的横坐标，直接以Hz为单位
xIndexFirst = linspace(0, fs, nFirst);
xIndexSecond = linspace(0, fs, nSecond);


harmIndex = getHarmIndex(fftWinDataFirst);  % 获取可能是谐波的数据角标，注意此时最小的角标是1，后续需要调整以符合定义
value = abs(fftWinDataFirst(harmIndex));  % 获取对应的幅值

% 获取基波幅值及对应角标
[fundaValue, fundaIndex] = max(sum(value, 1), [], 2);  % 就是找value里面最大的，返回列和最大的那一列的角标
% 删除基波对应的角标及幅值
RealfundaIndex = harmIndex(:, fundaIndex);  % 获取基波在fft变换结果里的角标
RealfundaValue = value(:, fundaIndex);  % 获取基波及其邻居对应的幅值
value(:, fundaIndex) = [];  % 删掉基波幅值
harmIndex(:, fundaIndex) = [];  % 删掉基波角标

% 画一下各次谐波的幅值，会发现有一些非常小，设一个阈值删去那些非常小的谐波（下一步完成）
xIndex = linspace(1, length(value), length(value));  % 横坐标为自然数
figure(8);
stem(xIndex, value(1, :));
hold on;
stem(xIndex, value(2, :), 'or');
hold on;
stem(xIndex, value(3, :));

% 删除幅值小于2的误判的谐波
delX = find(value(2,:)<2);
value(:, delX) = [];  % 删除
harmIndex(:, delX) = [];  % 删除

% 下面绘制删除误判谐波之后剩余的谐波，现在图里只剩下幅值不太小的数量很少的谐波幅值了
xIndex = linspace(1, length(value), length(value));  % 横坐标
figure(9);
stem(xIndex, value(1, :));
hold on;
stem(xIndex, value(2, :), 'or');
hold on;
stem(xIndex, value(3, :));

% 不做任何修正，计算谐波频率（就取幅值最大的，完全忽略左右邻居）
%format short g;
harmFreq = harmIndex./nFirst.*fs;
harmFreq(2, :);


% 修正谐波频率及幅值
[delta, vDelta] = blackmanHarrisFix(value);  % 计算修正量
fixHarmFreq = (harmIndex(2,:)+delta - 1)./nFirst.*fs;  % 这里要减1，因为fft结果角标应该从0开始
fixValue = (value(1,:)+2.*value(2,:)+value(3,:)).*vDelta./nFirst;

% 修正基波频率及幅值
[fDelta, fVDelta] = blackmanHarrisFix(RealfundaValue);
fixFundaFreq = (RealfundaIndex(2, :)+fDelta - 1)./nFirst.*fs;  % 这里要减一，因为fft结果角标应该从0开始
fixFundaValue = (RealfundaValue(1,:)+2.*RealfundaValue(2,:)+RealfundaValue(3,:)).*fVDelta./nFirst;

format longE;
HarmResult = [fixHarmFreq; fixValue]
FundaResult = [fixFundaFreq; fixFundaValue]




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