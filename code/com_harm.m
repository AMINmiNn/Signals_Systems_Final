fs = 10000;
data = csvread('1_1.csv');
%imf = emd(data(:,2),'Display',0);
%hht(imf,fs);
%根据希尔伯特-黄变换，得到跳变的具体时刻：0.4800。
%之后分段得到整体的基波幅值、频率、相角
[f,rms,phase] = myCal_FreFundamental(data(1:4801,2),fs)
[f,rms,phase] = myCal_FreFundamental(data(4802:9600,2),fs)
%原信号减去基波
t = (0:1/fs:0.4800)';
data(1:4801,2) = data(1:4801,2) - sqrt(2)*63.001363154034898*cos(50.998211896360928*2*pi*t-4.187023052672203);
t = (0.4801:1/fs:0.9599)';
data(4802:9600,2) = data(4802:9600,2) - sqrt(2)*75.603054533780139*cos(50.997671175724648*2*pi*(t-0.4801)-1.138025633128577);
%观察时域与频域图像
plot(data(:,2))
pspectrum(data(:,2),fs)
%发现4801点很异常，因此考虑取1~4800与4802~9600两段计算谐波，取平均值

%根据频谱构建对应谐波的带通滤波器
%1500Hz谐波
y_1500_1 = filter(filter_1500,data(1:4800,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
y_1500_2 = filter(filter_1500,data(4802:end,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
[f_1500_1,rms_1500_1] = myCal_FreFundamental(y_1500_1,fs)
[f_1500_2,rms_1500_2] = myCal_FreFundamental(y_1500_2,fs)
f_1500 = mean([f_1500_1,f_1500_2])
rms_1500 = mean([rms_1500_1,rms_1500_2])    

%1000Hz谐波
y_1000_1 = filter(filter_1000,data(1:4800,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
y_1000_2 = filter(filter_1000,data(4802:end,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
[f_1000_1,rms_1000_1] = myCal_FreFundamental(y_1000_1,fs)
[f_1000_2,rms_1000_2] = myCal_FreFundamental(y_1000_2,fs)
f_1000 = mean([f_1000_1,f_1000_2])
rms_1000 = mean([rms_1000_1,rms_1000_2])    

%650Hz
y_650_1 = filter(filter_650,data(1:4800,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
y_650_2 = filter(filter_650,data(4802:end,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
[f_650_1,rms_650_1] = myCal_FreFundamental(y_650_1,fs)
[f_650_2,rms_650_2] = myCal_FreFundamental(y_650_2,fs)
f_650 = mean([f_650_1,f_650_2])
rms_650 = mean([rms_650_1,rms_650_2])    

%550Hz
y_550_1 = filter(filter_550,data(1:4800,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
y_550_2 = filter(filter_550,data(4802:end,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
[f_550_1,rms_550_1] = myCal_FreFundamental(y_550_1,fs)
[f_550_2,rms_250_2] = myCal_FreFundamental(y_550_2,fs)
f_550 = mean([f_550_1,f_550_2])
rms_550 = mean([rms_550_1,rms_250_2])    


%250Hz
y_250_1 = filter(filter_250,data(1:4800,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
y_250_2 = filter(filter_250,data(4802:end,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
[f_250_1,rms_250_1] = myCal_FreFundamental(y_250_1,fs)
[f_250_2,rms_250_2] = myCal_FreFundamental(y_250_2,fs)
f_250 = mean([f_250_1,f_250_2])
rms_250 = mean([rms_250_1,rms_250_2])   

%250Hz
y_250_1 = filter(filter_250,data(1:4800,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
y_250_2 = filter(filter_250,data(4802:end,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
[f_250_1,rms_250_1] = myCal_FreFundamental(y_250_1,fs)
[f_250_2,rms_250_2] = myCal_FreFundamental(y_250_2,fs)
f_250 = mean([f_250_1,f_250_2])
rms_250 = mean([rms_250_1,rms_250_2])    

%230Hz
y_230_1 = filter(filter_230,data(1:4800,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
y_230_2 = filter(filter_230,data(4802:end,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
[f_230_1,rms_230_1] = myCal_FreFundamental(y_230_1,fs)
[f_230_2,rms_230_2] = myCal_FreFundamental(y_230_2,fs)
f_230 = mean([f_230_1,f_230_2])
rms_230 = mean([rms_230_1,rms_230_2])   


%200Hz
y_200_1 = filter(filter_200,data(1:4800,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
y_200_2 = filter(filter_200,data(4802:end,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
[f_200_1,rms_200_1] = myCal_FreFundamental(y_200_1,fs)
[f_200_2,rms_200_2] = myCal_FreFundamental(y_200_2,fs)
f_200 = mean([f_200_1,f_200_2])
rms_200 = mean([rms_200_1,rms_200_2])   

%150Hz
y_150_1 = filter(filter_150,data(1:4800,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
y_150_2 = filter(filter_150,data(4802:end,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
[f_150_1,rms_150_1] = myCal_FreFundamental(y_150_1,fs)
[f_150_2,rms_150_2] = myCal_FreFundamental(y_150_2,fs)
f_150 = mean([f_150_1,f_150_2])
rms_150 = mean([rms_150_1,rms_150_2])  

%100Hz
y_100_1 = filter(filter_100,data(1:4800,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
y_100_2 = filter(filter_100,data(4802:end,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
[f_100_1,rms_100_1] = myCal_FreFundamental(y_100_1,fs)
[f_100_2,rms_100_2] = myCal_FreFundamental(y_100_2,fs)
f_100 = mean([f_100_1,f_100_2])
rms_100 = mean([rms_100_1,rms_100_2])  

%70Hz
y_70_1 = filter(filter_70,data(1:4800,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
y_70_2 = filter(filter_70,data(4802:end,2));  % 直接使用设计好的滤波器进行滤波，filter函数是滤波函数
[f_70_1,rms_70_1] = myCal_FreFundamental(y_70_1,fs)
[f_70_2,rms_70_2] = myCal_FreFundamental(y_70_2,fs)
f_70 = mean([f_70_1,f_70_2])
rms_70 = mean([rms_70_1,rms_70_2])  