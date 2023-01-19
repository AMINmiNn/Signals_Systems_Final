%u基于Harr小波对必做信号的分析，求出了基波幅值随时间的变化
clc;
clear all;
wave = csvread('1_1.csv'); %load 数据
s1=wave(1:2048,2);
s2=wave(2049:4096,2);
s3=wave(4097:6144,2);
fs=10000; %采样率
r=1;
M=6;
n=2048;
syms phi;
phi=phi*2^(-M/2);
H=zeros(n,n*2^(-M));
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
f1=49.5;
w=2*pi*f1;
%G=[sin(w)*H,cos(w)*H];
G=zeros(n,2*n*2^(-M));
G1=zeros(n,n*2^(-M));
G2=zeros(n,n*2^(-M));
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

B1=inv(G'*G)*G'*s1;
a1=B1(1:n*2^(-M));
b1=B1(n*2^(-M)+1:2*n*2^(-M));
P1=H*a1;
Q1=H*b1;
A1=zeros(n,1);
for i=1:n
    A1(i,1)=(P1(i,1)^2+Q1(i,1)^2)^(1/2);
end


B2=inv(G'*G)*G'*s2;
a2=B2(1:n*2^(-M));
b2=B2(n*2^(-M)+1:2*n*2^(-M));
P2=H*a2;
Q2=H*b2;
A2=zeros(n,1);
for i=1:n
    A2(i,1)=(P2(i,1)^2+Q2(i,1)^2)^(1/2);
end

B3=inv(G'*G)*G'*s3;
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
A=cat(1,A1,A2,A3);
plot(t,A);