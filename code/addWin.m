clear all;
clc;

dataT = readmatrix('1_1.csv');
dataT = dataT.';
% 显示时域图
figure(1);
plot(dataT(1,:), dataT(2,:));
title('时域图');


fftResult = fft(dataT(2, :));
figure(2);
title('直接fft，对数幅值谱')
stem(dataT(1, :), 20*log10(abs(fftResult)));
figure(3);
stem(dataT(1, :), angle(fftResult));
title('直接fft相位谱')


Win = blackmanharris(length(dataT), 'periodic');
WinData = dataT(2, :).*Win';
figure(4);
plot(dataT(1, :), WinData(:));
title('加blackman窗，时域图')
fftWin = fft(WinData);
figure(5);
stem(dataT(1, :), 20*log10(abs(fftWin)));
title('加blackman窗，fft，对数幅值谱');
figure(6);
stem(dataT(1, :), angle(fftWin));
title('加blackman窗，fft，相位谱');



