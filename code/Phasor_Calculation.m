%��������

function [Phasor_funda,Phasor_harmo]=Phasor_Calculation(data)
% data��ԭʼ�����ź�,2�о���,��һ����ʱ��,�ڶ����ǵ�ѹ˲ʱֵ
%Phasor_funda�ǻ�����������������һ��3�еľ��󣬵�һ��ʱ�䣬�ڶ��з�ֵ����������λ
%Phasor_harmo��г������������һ��2�о��󣬵�һ����Ƶ�ʣ��ڶ����Ƿ�ֵ
t = data(:,1);
v = data(:,2);

%����Ƶ��
Fs = 10e3;

%ƫ�����ĵ�����ֵ
r0 = 0.5;  

%һ���ڵĲ�������
N0 = 200;   

%�ҷ�ֵ��Ծʱ��
t0 =  amptitude_step_find(v,N0,Fs); 

%г������
[Phasor_harmo, Funda] = anaHarmAllFromData(data,t0,Fs);

%�м����
FundaFrequency = Funda(1);
amptitude = zeros(2,1);
phase = zeros(2,1);

%�������
FundaAmp = zeros(2,1);
FundaPhase = zeros(2,1);

%ȷ����������λ
N0 = round(Fs/FundaFrequency);
M = floor(t0*Fs/N0);
i = 1;
j = 1;

%ͨ��ѭ����ȷ��й¶ʱ�ĳ������ȣ������ҵ��ȷ�Ļ�������λ
for N = N0*M:t0*Fs
    Xv = fft(v(1:N).*blackmanharris(N));
    [y2,index2] = max(abs(Xv));  %����ֵ����λ��
    y1 = abs(Xv(index2 + 1));
    alpha = y2/y1;
    myfun = @(r) deviation(r,alpha);
    r = fzero(myfun,r0);  %ƫ����r
    amptitude(j)= 2*y2*pi*r*(1-r^2)*(4-r^2)*(9-r^2)/(sin(r*pi)*(12.915-1.22511*r^2 ...
        +0.02913*r^4-0.00006*r^6))/(N)/sqrt(2);
    phase(j) = angle(Xv(index2)) - r*pi;
    j = j+1;
end
    [~,index0] = max(amptitude);
    FundaPhase0 = phase(index0);

%�Ӵ���ÿ����ķ�ֵ
L = N0*4;  %�趨���Ŀ��
for  j = L/2+1:t(end)*Fs-L/2
     v1 = v(j-L/2:j+L/2-1);  %���Ӵ����ź�
     v2 = v1.*blackmanharris(L);   %�Ӵ�����֮��
     Xv2 = fft(v2);
    [y2,index2] = max(abs(Xv2));  %����ֵ����λ��
    y1 = abs(Xv2(index2 + 1));
    alpha = y2/y1;
    myfun = @(r) deviation(r,alpha);
    r = fzero(myfun,r0);  %ƫ����r    
    FundaAmp(j) = 2*y2*pi*r*(1-r^2)*(4-r^2)*(9-r^2)/(sin(r*pi)*(12.915-1.22511*r^2 ...
        +0.02913*r^4-0.00006*r^6))/L;
end

%��ֵ����2�����������Ȱ�����Ծǰ���ְ�����Ծ��Ĳ��֣���Ծ��ǰ�ķ�ֵȡAmp2(t0*Fs-L/2)����Ծ�Ժ�ķ�ֵȡAmp2(t0*+L/2)��
FundaAmp(t0*Fs-L/2+1:t0*Fs-1) = FundaAmp(t0*Fs-L/2);
FundaAmp(t0*Fs:t0*Fs+L/2) = FundaAmp(t0*Fs+L/2);
FundaAmp1 = FundaAmp/sqrt(2);


%��λ����2��������Ծ��������λ�����Լ�ȡ�����λ
FundaPhase = FundaPhase0 + 2*pi*(FundaFrequency-50)*t(1:(end-L/2-1));
FundaPhase1 = mod(FundaPhase,2*pi);
FundaPhase1 = FundaPhase1.*(0<=FundaPhase1 & FundaPhase1 <= pi) + (FundaPhase1 - 2*pi).*(pi<FundaPhase1 & FundaPhase1<2*2*pi);  

figure(1)
subplot(3,1,1)
amp = FundaAmp.*cos(FundaPhase + 2*pi*50*t(1:(end-L/2-1)));
plot(t(L/2+1:length(amp)),amp(L/2+1:end));
hold on;
plot(t(L/2+1:length(amp)),v(L/2+1:length(amp)));
xlabel("ʱ��/s");
ylabel("��ֵ/V");
legend("����ֵ","ʵ��ֵ");
title("����������ֵ����λ");
hold off;

subplot(3,1,2)
plot(t(L/2+1:end-L/2-1),FundaAmp1(L/2+1:end));
title("������ֵ����");
xlabel("ʱ��/s");
ylabel("��Чֵ/V");

subplot(3,1,3)
plot(t(L/2+1:end-L/2-1),FundaPhase1(L/2+1:end));
title("������λ����");
xlabel("ʱ��/s");
ylabel("�����λ");

Phasor_funda = [t(L/2+1:end-L/2-1),FundaAmp1(L/2+1:end),FundaPhase1(L/2+1:end)];

    function f = deviation(r,alpha)
    f = alpha*(r+3)*(2*r^6-12*r^5-941*r^4+3844*r^3+35041*r^2-77802*r-390632)+...
        (2*r^6-971*r^4+40837*r^2-430500)*(r-4);
    end

    function t0 = amptitude_step_find(v,N0,Fs)
        imf = emd(v);%���������ģ̬����
        [~, ~, t1, imfinsf] = hht(imf(:,1),Fs);
        [~, index] = max(imfinsf(2*N0:end-2*N0));  %Ĭ��ǰ�����ܲ������ֽ�Ծ
        t0 = t1(index + 2*N0 - 1); %ͨ��ϣ�����ر任��ȡ��Ծʱ��
    end

    
    function [Harm, Funda] = anaHarmAllFromData(data,t0,Fs)
    dataT = data.';
    fs = 10e3;  % ����Ƶ��
    [peakIndex, part1, part2] = getPeakIndex(dataT,t0,Fs);  % ��ֳ������֣�ֻȡ������
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
    

    % ����ԭʼ��fft
    fftDataFirst = fft(dataFirst(2, :));
    fftDataSecond = fft(dataSecond(2, :));

    % ��Ծǰ�ķ������
    [HarmResult1, FundaResult1] = anaHarm(fftWinDataFirst, fftDataFirst);
    % ��Ծ��ķ������
    %[HarmResult2, FundaResult2] = anaHarm(fftWinDataSecond, fftDataSecond);
    % ����Ծǰ��Ƶ�ʣ���ֵд���ļ�
   
    Harm = HarmResult1(2:3, :);
    Harm = [Harm(1, :); Harm(2, :)./sqrt(2)];
    Funda = FundaResult1(2:3, :);
    Funda = [Funda(1, :); Funda(2, :)./sqrt(2)];
    end


    function [HarmResult, FundaResult] = anaHarm(fftWinData, fftData)
        fs = 10e3;  % ����Ƶ��
        dataN = length(fftWinData);  % ��ȡ����
        harmIndex = getHarmIndex(fftWinData);  % ��ȡг���Ǳ꣨��ʱ��г����Ȼ����������
        value = abs(fftWinData(harmIndex));  % ��ȡг����ֵ
        [fundaValue, fundaIndex] = max(sum(value, 1), [], 2);  % ȡ�������Ǳ�
        RealfundaIndex = harmIndex(:, fundaIndex);  % ��ȡ������fft�еĽǱ�
        RealfundaValue = value(:, fundaIndex);  % ��ȡ������ֵ
        value(:, fundaIndex) = [];  % ɾ��������ֵ
        harmIndex(:, fundaIndex) = [];  % ɾ�������Ǳ�
        
        % ��ͼ����г���������ھӶ����������ᷢ�������еģ��еķ�����ֵ��С
        xIndex = linspace(1, length(value), length(value));
        figure(2);
        stem(xIndex, value(1, :));
        hold on;
        stem(xIndex, value(2, :), 'or');
        hold on;
        stem(xIndex, value(3, :));
        title('δ��������ǰ��г��')
        
        % �趨��ֵΪ2������ֵС��2��һ����Ϊ�����У�ɾȥ���е�г��
        level = 2;
        delX = find(value(2,:)<level);
        value(:, delX) = [];  % ɾ��
        harmIndex(:, delX) = [];  % ɾ�� 
        
        % ���»�ͼ��������г�����������ھӣ���ʱû�����еģ�����ʣ���г����������
        xIndex = linspace(1, length(value), length(value));  % ������
        figure(3);
        stem(xIndex, value(1, :));
        hold on;
        stem(xIndex, value(2, :), 'or');
        hold on;
        stem(xIndex, value(3, :));
        title('�������к��г��')
        
        % ����г��Ƶ�ʼ���ֵ������blackmanharris�����������̼���
        [delta, vDelta] = blackmanHarrisFix(value);  % ����������
        fixHarmFreq = (harmIndex(2,:)+delta - 1)./dataN.*fs;  % ����Ҫ��1����Ϊfft����Ǳ�Ӧ�ô�0��ʼ
        fixValue = (value(1,:)+2.*value(2,:)+value(3,:)).*vDelta./dataN;
    
        % ��������Ƶ�ʼ���ֵ
        [fDelta, fVDelta] = blackmanHarrisFix(RealfundaValue);
        fixFundaFreq = (RealfundaIndex(2, :)+fDelta - 1)./dataN.*fs;  % ����Ҫ��һ����Ϊfft����Ǳ�Ӧ�ô�0��ʼ
        fixFundaValue = (RealfundaValue(1,:)+2.*RealfundaValue(2,:)+RealfundaValue(3,:)).*fVDelta./dataN;
        
        % ��ȡԭʼ��λ
        rawAngle = angle(fftData);
        % �����������г����λ
        harmAngle = rawAngle(harmIndex(2, :))  - delta.*pi;
        % ����������Ļ�����λ
        fundaAngle = rawAngle(RealfundaIndex(2, :))  - fDelta.*pi;
    
    
        format longE;
        HarmResult = [harmIndex(2, :); fixHarmFreq; fixValue; harmAngle];  % ԭʼ�Ǳ� ����Ƶ�� ������ֵ ������λ
        FundaResult = [RealfundaIndex(2, :); fixFundaFreq; fixFundaValue; fundaAngle];
    end
    
    
    % blackmanHarris�����β�ֵ��������
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
    function [peakIndex, part1, part2] = getPeakIndex(dataT,t0,Fs)
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
            if peakIndex(i) >= t0*Fs
                cutIndex = i;
                break;
            end
        end
        part1 = peakIndex(1:cutIndex-1);
        part2 = peakIndex(cutIndex:end);
    end

    writematrix([Phasor_funda,[Phasor_harmo';zeros(size(Phasor_funda,1)-size(Phasor_harmo,2),2)]],'1_1������֤.csv');
end
