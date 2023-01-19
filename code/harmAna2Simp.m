clf;

% �ֳɽ�Ծǰ�ͽ�Ծ�������ֽ��з���
dataT = readmatrix('1_1.csv');
dataT = dataT.';
fs = 10e3;  % ����Ƶ��
[peakIndex, part1, part2] = getPeakIndex(dataT);  % ��ֳ������֣�ֻȡ������
dataFirst = dataT(:, part1(1):part1(end));
dataSecond = dataT(:, part2(1):part2(end));


% �������г���
nFirst = length(dataFirst);  
nSecond = length(dataSecond);
% ����blackmanharris��
winFirst = blackmanharris(nFirst, 'periodic');
winSecond = blackmanharris(nSecond, 'periodic');
% ����Ӵ����ʱ������
winDataFirst = dataFirst(2, :).*winFirst';
winDataSecond = dataSecond(2, :).*winSecond';
% �ԼӴ���Ľ������fft
fftWinDataFirst = fft(winDataFirst);
fftWinDataSecond = fft(winDataSecond);

% ����fft�ĺ����ֱ꣬����HzΪ��λ
xIndexFirst = linspace(0, fs, nFirst);
xIndexSecond = linspace(0, fs, nSecond);


harmIndex = getHarmIndex(fftWinDataFirst);  % ��ȡ������г�������ݽǱ꣬ע���ʱ��С�ĽǱ���1��������Ҫ�����Է��϶���
value = abs(fftWinDataFirst(harmIndex));  % ��ȡ��Ӧ�ķ�ֵ

% ��ȡ������ֵ����Ӧ�Ǳ�
[fundaValue, fundaIndex] = max(sum(value, 1), [], 2);  % ������value�������ģ������к�������һ�еĽǱ�
% ɾ��������Ӧ�ĽǱ꼰��ֵ
RealfundaIndex = harmIndex(:, fundaIndex);  % ��ȡ������fft�任�����ĽǱ�
RealfundaValue = value(:, fundaIndex);  % ��ȡ���������ھӶ�Ӧ�ķ�ֵ
value(:, fundaIndex) = [];  % ɾ��������ֵ
harmIndex(:, fundaIndex) = [];  % ɾ�������Ǳ�

% ��һ�¸���г���ķ�ֵ���ᷢ����һЩ�ǳ�С����һ����ֵɾȥ��Щ�ǳ�С��г������һ����ɣ�
xIndex = linspace(1, length(value), length(value));  % ������Ϊ��Ȼ��
figure(8);
stem(xIndex, value(1, :));
hold on;
stem(xIndex, value(2, :), 'or');
hold on;
stem(xIndex, value(3, :));

% ɾ����ֵС��2�����е�г��
delX = find(value(2,:)<2);
value(:, delX) = [];  % ɾ��
harmIndex(:, delX) = [];  % ɾ��

% �������ɾ������г��֮��ʣ���г��������ͼ��ֻʣ�·�ֵ��̫С���������ٵ�г����ֵ��
xIndex = linspace(1, length(value), length(value));  % ������
figure(9);
stem(xIndex, value(1, :));
hold on;
stem(xIndex, value(2, :), 'or');
hold on;
stem(xIndex, value(3, :));

% �����κ�����������г��Ƶ�ʣ���ȡ��ֵ���ģ���ȫ���������ھӣ�
%format short g;
harmFreq = harmIndex./nFirst.*fs;
harmFreq(2, :);


% ����г��Ƶ�ʼ���ֵ
[delta, vDelta] = blackmanHarrisFix(value);  % ����������
fixHarmFreq = (harmIndex(2,:)+delta - 1)./nFirst.*fs;  % ����Ҫ��1����Ϊfft����Ǳ�Ӧ�ô�0��ʼ
fixValue = (value(1,:)+2.*value(2,:)+value(3,:)).*vDelta./nFirst;

% ��������Ƶ�ʼ���ֵ
[fDelta, fVDelta] = blackmanHarrisFix(RealfundaValue);
fixFundaFreq = (RealfundaIndex(2, :)+fDelta - 1)./nFirst.*fs;  % ����Ҫ��һ����Ϊfft����Ǳ�Ӧ�ô�0��ʼ
fixFundaValue = (RealfundaValue(1,:)+2.*RealfundaValue(2,:)+RealfundaValue(3,:)).*fVDelta./nFirst;

format longE;
HarmResult = [fixHarmFreq; fixValue]
FundaResult = [fixFundaFreq; fixFundaValue]




function [delta, vDelta] = blackmanHarrisFix(y)
    alpha = (y(3,:)-y(1,:))./y(2,:);
    delta = 0.93891885.*alpha - 0.08203810.*alpha.^3 + 0.01541147.*alpha.^5 - 0.00316988.*alpha.^7;
    vDelta = 1.65866360 + 0.44865377.*delta.^2 + 0.06478728.*delta.^4 + 0.00693833.*delta.^6;
end






% ������ֵ�������������δ��ֵ�ĽǱ�
function harmIndex = getHarmIndex(fftData)
    fftData = abs(fftData);
    len = length(fftData)/2;  % ��ȡ���ȣ�ȡһ�뼴�ɣ���Ϊ�ǶԳƵ�
    harmIndex = zeros(3, 5);  % ��һ��Ϊ���ھӣ�������Ϊ���ھ�
    index = 1;  % �洢�õĽǱ�
    for i = 1 + 1: len - 1
        if fftData(i) > fftData(i-1) && fftData(i) > fftData(i+1)
            harmIndex(1, index) = i-1;
            harmIndex(2, index) = i;
            harmIndex(3, index) = i+1;
            index = index + 1;  % �Ǳ�����
        end
    end
end







% �ҵ����ҷ�ֵ��Ӧ�ĽǱ�
function [peakIndex, part1, part2] = getPeakIndex(dataT)
    %peakIndex ���������з�ֵ�ĽǱ꣬part1����������δ��Ծ�ķ�ֵ�ĽǱ꣬part2���������н�Ծ��ķ�ֵ�ĽǱ�
    len = length(dataT(1, :));  % ��ȡʱ�򳤶�
    peakIndex = zeros(1, 10);  % �洢��ֵ�Ǳ�
    saveIndex = 1;  % �����洢�ĽǱ�
    margin = 10;  % �۲��򣬷�ֵҪ�����Ҹ�margin�����
    for i = 1+margin: len-margin
        if dataT(2, i) > max(dataT(2, i-margin: i-1)) && dataT(2, i) > max(dataT(2, i+1:i+margin))
            peakIndex(saveIndex) = i;
            saveIndex = saveIndex + 1;   % �Ǳ�����1
        end
    end
    
    cutIndex = 0;  % �ָ���Ծǰ�ͽ�Ծ��ĽǱ�
    for i = 2:length(peakIndex)  % ��ǰһ����ֵ��1.1��С�ں�һ����ֵ����������ֵ������ͬ���ڴ˴����л���
        if dataT(2, peakIndex(i)) > dataT(2, peakIndex(i-1)).*1.1
            cutIndex = i;
            break;
        end
    end
    part1 = peakIndex(1:cutIndex-1);
    part2 = peakIndex(cutIndex:end);
end