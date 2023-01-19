M=csvread('1_1.csv');
fs=10000;
len=length(M);
window=600;  % 可以修改 window值 观察效果
k_range=1:(len-window);
f=zeros(len,1);
A=zeros(len,1);
phase=zeros(len,1);
for k=k_range
    [f(k),A(k),phase(k)]=myCal_FreFundamental(M(k:k+window,2),M(k:k+window,1),fs);
    phase(k)=mod(phase(k)-50*2*pi*M(k,1),2*pi);
end
%% 该部分代码是为补全 A和phase 并不必要
A(len-window+1:end)=A(len-window);
for j =len-window+1:len
    phase(j)=mod(phase(len-window)+2*pi*f(len-window)*(M(j,1)-M(len-window,1))-2*pi*50*(M(j,1)-M(len-window,1)),2*pi);
end

%% 绘图
subplot(2,1,1)
plot(M(1:len-window,1),A(1:len-window));
title('RMS')
hold on;
subplot(2,1,2)
plot(M(1:len-window,1),phase(1:len-window));
title('phase')
%% write csv
%%writematrix([M(:,1) A phase],'1_complusory_fundamental.csv')

