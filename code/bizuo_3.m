clc;close all;
%小波变换计算阶跃时间

table = csvread('1_1.csv'); %load 数据
V1=table(1:2048,2);
V2=table(2049:4096,2);
V3=table(4097:6144,2);
fs=10000; %采样�?
r=1;
M=6;
n=2048;
syms phi;
phi=phi*2^(-M/2);
H=zeros(n,n*2^(-M));%2048*32矩阵
for i=1:n
    for j=1:n*2^(-M)
        if (2^(-M)*(i-1)-(j-1)>=0 && 2^(-M)*(i-1)-(j-1)<1)
            phi=1;
        else
            phi=0;
        end
        H(i,j)=phi*2^(-M/2);
    end
end
f1=49;
w=2*pi*f1;
%G=[sin(w)*H,cos(w)*H];
G=zeros(n,2*n*2^(-M));%2048*64矩阵，容纳G1和G2
G1=zeros(n,n*2^(-M));%2048*32矩阵
G2=zeros(n,n*2^(-M));%2048*32矩阵
for i=1:n
    G1(i,:)=H(i,:)*sin(w*(i-1)/fs);
    G2(i,:)=H(i,:)*cos(w*(i-1)/fs);
end
for j=1:2*n*2^(-M)
    if j<=32
        G(:,j)=G1(:,j);
    else
        G(:,j)=G2(:,j-32);
    end
end

B1=inv(G'*G)*G'*V1;
a1=B1(1:n*2^(-M));
b1=B1(n*2^(-M)+1:2*n*2^(-M));
P1=H*a1;
Q1=H*b1;
A1=zeros(n,1);
for i=1:n
    A1(i,1)=(P1(i,1)^2+Q1(i,1)^2)^(1/2);
end

B2=inv(G'*G)*G'*V2;
a2=B2(1:n*2^(-M));
b2=B2(n*2^(-M)+1:2*n*2^(-M));
P2=H*a2;
Q2=H*b2;
A2=zeros(n,1);
for i=1:n
    A2(i,1)=(P2(i,1)^2+Q2(i,1)^2)^(1/2);
end

B3=inv(G'*G)*G'*V3;
a3=B3(1:n*2^(-M));
b3=B3(n*2^(-M)+1:2*n*2^(-M));
P3=H*a3;
Q3=H*b3;
A3=zeros(n,1);
for i=1:n
    A3(i,1)=(P3(i,1)^2+Q3(i,1)^2)^(1/2);
end

T=0:6143;
t=T/fs;
A=cat(1,A1,A2,A3);%按列连接小波变换幅�??
plot(t,A);

%分离谐波计算阶跃时间
table = xlsread('1_1.csv');
time = table(1:end,1);
val = table(1:end,2);

s=table(1:9600,2);
fs=10000; 
ts = 0:0.0001:0.9599;
s1=87.6812*sin(2*pi*49*ts+210/180*pi)+0.0510*sin(2*pi*97.4166*ts-98.3307/180*pi)...
+0.0484*sin(2*pi*247*ts+108.1096/180*pi)+0.00003*sin(2*pi*292.0768*ts+91.8465/180*pi)...
+0.0139*sin(2*pi*340.4348*ts-27.7664/180*pi)+0.01580*sin(2*pi*445*ts+326.443/180*pi);

s2=s-s1';
N=9600;
n=0:N-1;
t=n/fs;
figure(4);
plot(t,s2);