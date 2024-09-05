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
% [Lh1,Lh2,Lh3]=Getlambda1(lambdav,posi,c1_w,c2_w,c3_w,n,M);
% Idx=[3 3 1 2 1 1 1];
% lambda1=lambdav;
% [row col]=size(lambda1);
% Ck=n;
% posv=zeros(Ck,1);
% lambda2=zeros(2,Ck);
% lh1=zeros(2,col); i1=1;rat11=zeros(2,col);
% lh2=zeros(2,col); i2=1;rat22=zeros(2,col);
% lh3=zeros(2,col); i3=1;rat33=zeros(2,col);
% for ix=1:Ck
%     v=[];
%     Lhv0=0;
%     Lhv1=0;
%     lambda2(:,ix)=0;
%     for iy=1:col 
%           if(Idx(iy)==ix)
%              v=[v iy];
%              a_j0=lambda1(1,iy);
%              a_j1=lambda1(2,iy);
%              Lhv00=(1/(sqrt(2*pi)*1)*exp(-1/(2*1)*abs(a_j0*d1(locs(iy))-d2(locs(iy)))^2/(1+a_j0^2)));
%              Lhv11=(1/(sqrt(2*pi)*1)*exp(-1/(2*1)*abs(a_j1*d1(locs(iy))-d3(locs(iy)))^2/(1+a_j1^2)));
%              if(Lhv00>Lhv0)
%                 lambda2(1,ix) = lambda1(1,iy);
%                 Lhv0=Lhv00;
%              end
%              if(Lhv11>Lhv1)
%                 lambda2(2,ix) = lambda1(2,iy);
%                 Lhv1=Lhv11;
%              end
%              if(ix==1)
%                  lh1(1,i1)= Lhv00;
%                  lh1(2,i1)= Lhv11;
%                  rat11(1,i1)=a_j0;
%                  rat11(2,i1)=a_j1;
%                  i1=i1+1;
%              end
%              if(ix==2)
%                  lh2(1,i2)= Lhv00;
%                  lh2(2,i2)= Lhv11;
%                  rat22(1,i2)=a_j0;
%                  rat22(2,i2)=a_j1;
%                  i2=i2+1;
%              end
%              if(ix==3)
%                  lh3(1,i3)= Lhv00;
%                  lh3(2,i3)= Lhv11;
%                  rat33(1,i3)=a_j0;
%                  rat33(2,i3)=a_j1;
%                  i3=i3+1;
%              end             
%           end
%     end
%     posv(ix,1) = min(v);
% end
% lambda = zeros(2,Ck);
% sposv=sort(posv);
% posv1=zeros(Ck,1);
% for ix=1:Ck;
%     posv1(ix)=find(sposv== posv(ix));
% end
% for ix=1:Ck;
%     iy=find(posv1==ix);
%     lambda(:,ix)=lambda2(:,iy);
% end
%=========================================================================%
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

