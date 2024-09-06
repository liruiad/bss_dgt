clc;
clear all;
close();
L=1024;
Ls=1024*15;%25600*2;
fs=50000;
m=3;n=3;
N=L;
N0=1;
load s1 s1;
load s2 s2;
load s3 s3;
load c1_w c1_w;
load c2_w c2_w;
load c3_w c3_w;
load x1 x1;
load x2 x2;
load x3 x3;
load X  X;
load A A;

M=(Ls+L)/N0;
subplot(131);
contour(abs(c1_w(:,1:N/2)'));
subplot(132);
contour(abs(c2_w(:,1:N/2)'));
subplot(133);
contour(abs(c3_w(:,1:N/2)'));
%%=====================================================================%%
d1=zeros(N,1);
d2=zeros(N,1);
d3=zeros(N,1); 
var1=zeros(N,1);
var2=zeros(N,1);
var3=zeros(N,1);
for ix=1:N/2;
    d1(ix)=sum(c1_w(:,ix))/M;
    d2(ix)=sum(c2_w(:,ix))/M;
    d3(ix)=sum(c3_w(:,ix))/M;
    var1(ix)=var(c1_w(:,ix));
    var2(ix)=var(c2_w(:,ix));
    var3(ix)=var(c3_w(:,ix));
end

val = abs(d1(1:N/2))+abs(d2(1:N/2))+abs(d3(1:N/2));
val(1)=0;
val=val.^4;
val=(mapminmax(val',0,1))';
[pks11,locs123] = findpeaks(val,'minPeakHeight',0.25,'MinPeakDistance',1);% 0.9 10
valv=sort(val,'descend');

lenp=length(locs123);
lambda0=zeros(2,lenp);
lambda1=[];
Lh=[];
% Lh1=[];
locs=[];
theta=1;
for ix=1:lenp 
        lambda0(1,ix)=abs(d2(locs123(ix))/d1(locs123(ix)));
        lambda0(2,ix)=abs(d3(locs123(ix))/d1(locs123(ix)));
        lambda1=[lambda1 lambda0(:,ix)];
        locs=[locs locs123(ix)];
        zz=abs(lambda0(1,ix)*d1(locs123(ix))-d2(locs123(ix)))^2 + abs(lambda0(2,ix)*d1(locs123(ix))-d3(locs123(ix)))^2 ...
           +abs(lambda0(2,ix)*d2(locs123(ix))-lambda0(1,ix)*d3(locs123(ix)))^2;
        Lh=[Lh 1/(sqrt(2*pi)^3*theta^3)*exp(-1/(2*theta^2)*zz/(1+abs(lambda0(1,ix))^2+abs(lambda0(2,ix))^2))];
end
posi=[];
Hvhv=[];
lambdav=[];
Lh=(Lh-min(Lh))/(max(Lh)-min(Lh));
for ix=1:length(Lh)
    if(Lh(ix)>0.8) %0.8
        posi=[posi locs(ix)];
        lambdav=[lambdav lambda1(:,ix)];
        Hvhv=[Hvhv Lh(ix)];
    end
end

[Idx,lambda,lh1,lh2,lh3,rat11,rat22,rat33]=Getlambda(lambdav,posi,d1,d2,d3,n);
A1=[1 1 1;lambda];
SR=inv(A1)*X;
save SR SR;
sr1=SR(3,:);
sr2=SR(2,:);
sr3=SR(1 ,:); 
% 
s1=(mapminmax(s1',-1,1))';
s2=(mapminmax(s2',-1,1))';
s3=(mapminmax(s3',-1,1))';
sr1=mapminmax(sr1,-1,1);
sr2=mapminmax(sr2,-1,1);
sr3=mapminmax(sr3,-1,1);
% -----------%
figure;
subplot(321)
plot(s1);
xlabel('s1');
subplot(322)
plot(sr1);
xlabel('sr1');
subplot(323)
plot(s2);
xlabel('s2');
subplot(324)
plot(sr2);
xlabel('sr2');
subplot(325)
plot(s3);
xlabel('s3');
subplot(326)
plot(sr3);
xlabel('sr3');
% 
figure;
subplot(311)
plot(s1);
subplot(312)
plot(s2);
subplot(313)
plot(s3);

figure;
subplot(311)
plot(sr1);
subplot(312)
plot(sr2);
subplot(313)
plot(sr3);

