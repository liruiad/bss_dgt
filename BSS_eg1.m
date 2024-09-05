clc;
clear all;close all;
L=1024;
g=zeros(L,1);
for ix=0:L-1;
    g(ix+1)=0.1*exp(-pi/2*((ix-(L-1)/2)*0.001).^2);
end
%%=======================================================================%%
Ls=1024*1;fs=1000;
load xBper_i xBper_i;
load xBper_o xBper_o;
%%=======================================================================%%
m=2;n=2;% the number of observed signals
%A=2*rand(m,n);
A= [1.1769,0.8878;0.4940,1.7366];
s1=xBper_i(1:Ls);  % source signal s1
s2=xBper_o(1:Ls);  % source signal s2
S=[s1;s2];
subplot(211)
plot(s1);
subplot(212)
plot(s2);


X=A*S;     % mixing signal 
x1=X(1,:); % mixing signal x1
x2=X(2,:); % mixing signal x2

figure;
subplot(211)
plot(x1);
subplot(212)
plot(x2);

N=L;N0=1;
c1_w=DGT(L,N,N0,x1,g);
c2_w=DGT(L,N,N0,x2,g);
figure;
subplot(121);
contour(abs(c1_w(:,1:N/2)'));
subplot(122);
contour(abs(c2_w(:,1:N/2)'));

M=(Ls+L)/N0;
d1=zeros(N,1);
d2=zeros(N,1);
for ix=1:N/2;
    d1(ix)=sum(c1_w(:,ix))/M; % eq.(20) eq(26)
    d2(ix)=sum(c2_w(:,ix))/M; % eq.(20) eq(26)
end
val = abs(d1(1:N/2))+abs(d2(1:N/2));
val(1)=0;
val=val.^2;
val=(mapminmax(val',0,1))';
[pks11,locs123] = findpeaks(val,'minPeakHeight',0.2,'MinPeakDistance',1);% 0.9 10
valv=sort(val);
figure;plot(val);
p_f_m=zeros(3,length(locs123));
p_f_m(1,:)=locs123';
p_f_m(2,:)=round(locs123'./1024*20000,1);
p_f_m(3,:)=pks11';

lenp=length(locs123);
lambda0=zeros(1,lenp);
lambda1=[];
Lh=[];
% Lh1=[];
locs=[];
theta=1;
for ix=1:lenp
        lambda0(1,ix)=abs(d2(locs123(ix))/d1(locs123(ix))); % eq.(21) eq(27)
        lambda1=[lambda1 lambda0(:,ix)];
        locs=[locs locs123(ix)];
        zz=abs(lambda0(1,ix)*d1(locs123(ix))-d2(locs123(ix)))^2;
        Lh=[Lh 1/(sqrt(2*pi)^2*theta^2)*exp(-1/(2*theta^2)*zz/(1+abs(lambda0(1,ix))^2))]; % eq.(34)
end
posi=[];
Hvhv=[];
lambdav=[];
Lh=(Lh-min(Lh))/(max(Lh)-min(Lh));
for ix=1:length(Lh)
    if(Lh(ix)>0.95) %0.8
        posi=[posi locs(ix)];
        lambdav=[lambdav lambda1(:,ix)];
        Hvhv=[Hvhv Lh(ix)];
    end
end
Lh_locs=[Lh;locs123'];
Hvhv_posi=[Hvhv;posi];
[Idx,C]=kmeans(lambdav',n);
%Idx= [1;1;1;1;1;1;1;1;2;2;2;2]; L=256;
%Idx= [2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;1;2;1;2;1;1;1;1;1]; L=1024;
lambda=zeros(1,2);
Lhv0=0;
Lhv1=0;
 for iy=1:length(lambdav)
    if(Idx(iy)==1)
        if(Hvhv(iy)>Lhv0)
            Lhv0 = Hvhv(iy);
            lambda(1,1)=lambdav(iy);
        end
    end
    if(Idx(iy)==2)
        if(Hvhv(iy)>Lhv1)
            Lhv1 = Hvhv(iy);
            lambda(1,2)=lambdav(iy);
        end
    end
end

A=[1 1;lambda];
SR=inv(A)*X;
sr1=SR(2,:)';
sr2=SR(1,:)';
% 
s11=mapminmax(s1,-1,1);
s21=mapminmax(s2,-1,1);
sr1=mapminmax(sr1',-1,1);
sr2=mapminmax(sr2',-1,1);
% -----------%
figure;
subplot(221)
plot(s11);
xlabel('s11');
subplot(222)
plot(sr1);
xlabel('sr1');
subplot(223)
plot(s21);
xlabel('s21');
subplot(224)
plot(sr2);
xlabel('sr2');


Y1=abs((fft(s1))).^2;
Y2=abs(fft(s2)).^2;
Y1=mapminmax(Y1,0,1);
Y2=mapminmax(Y2,0,1);
figure; 
subplot(211);
plot(Y1(1:Ls/2));
subplot(212);
plot(Y2(1:Ls/2));
[pks1,locs1] = findpeaks(Y1(1:Ls/2),'minPeakHeight',0.1,'MinPeakDistance',1);
[pks2,locs2] = findpeaks(Y2(1:Ls/2),'minPeakHeight',0.1,'MinPeakDistance',1);

figure;
subplot(211)
plot(sr1);
subplot(212)
plot(sr2);
