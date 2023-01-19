clear all;
clc;

dataT = readmatrix('1_1.csv');
dataT = dataT.';
% ��ʾʱ��ͼ
figure(1);
plot(dataT(1,:), dataT(2,:));
title('ʱ��ͼ');


fftResult = fft(dataT(2, :));
figure(2);
title('ֱ��fft��������ֵ��')
stem(dataT(1, :), 20*log10(abs(fftResult)));
figure(3);
stem(dataT(1, :), angle(fftResult));
title('ֱ��fft��λ��')


Win = blackmanharris(length(dataT), 'periodic');
WinData = dataT(2, :).*Win';
figure(4);
plot(dataT(1, :), WinData(:));
title('��blackman����ʱ��ͼ')
fftWin = fft(WinData);
figure(5);
stem(dataT(1, :), 20*log10(abs(fftWin)));
title('��blackman����fft��������ֵ��');
figure(6);
stem(dataT(1, :), angle(fftWin));
title('��blackman����fft����λ��');



