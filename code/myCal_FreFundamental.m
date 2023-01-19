function [f0,A,phase] = myCal_FreFundamental(x,fs)
    N = length(x);
    %% 获得hamming窗的特性曲线，alpha与dmu的关系
    % 建立alpha-dmu曲线
    % This curve demonstrates the relationship between error Dmu and alpha
    w = hamming(N);
    m = 500; 
    Dmu_li = 0:1/m:0.5;
    W = czt(w, m + 1, exp(-j * 2 * pi / (m * N)), 1);
    index = 1:length(Dmu_li);
    Alpha_li = abs(W(end - index + 1)) ./ abs(W(index));

    %% 计算加窗后的DFT:x*w
    X = fft(x.*w);
    m = [0:length(X) - 1] / length(X);

    %% 查找加窗DFT最大次大值
    Xa = abs(X);
    m0 = find(Xa == max(Xa), 1);
    if Xa(m0 + 1) > Xa(m0 - 1)
        m1 = m0 + 1;
    else
        m1 = m0 - 1;
    end

    %% 频率估计
    alpha = Xa(m1) / Xa(m0);
    dmu = interp1(Alpha_li, Dmu_li, alpha);

    if m1 < m0
        dmu = -dmu;
    end

    f1 = (m0 + dmu - 1) / N;    f0=f1 * fs;
    % disp(["fundamental frequency:",f0]);
    %% 幅度(RMS)和相位估计
    k = 0:N-1;
    W1 = sum(w'.*exp(j * 2 * pi * dmu * k / N));
    A = abs(2 * X(m0) / W1/sqrt(2));
    % disp(["amplitude", A]);
    phase = angle(X(m0)) - angle(W1);
    % disp(["phase", phase]);

    
end


