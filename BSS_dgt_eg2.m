 clc;
close();
clear all;
L=1024;
for ix=0:L-1;
    g(ix+1)=0.1*exp(-pi/2*((ix-(L-1)/2)*0.001).^2);
end
save g g;
%%=======================================================================%%
load('D:\Code\BSS_DGT\JNU-Bearing-Dataset-main\data_all');
Ls=1024*15;%25600*2;
s1=ob10002(1:Ls);
s2=ib10002(1:Ls);
s3=tb10002(1:Ls);
fs=50000;
save s1 s1;
save s2 s2;
save s3 s3;
subplot(311)
plot(s1);
subplot(312)
plot(s2);
subplot(313)
plot(s3);
%%=======================================================================%%
m=3;n=3;% the number of observed signals
% A=rand(m,n); 
A= [0.430207391329584,0.979748378356085,0.258064695912067;
 0.184816320124136,0.438869973126103,0.408719846112552;
 0.904880968679893,0.111119223440599,0.594896074008614];
S=[s1';s2';s3'];
X=A*S;

x1=X(1,:);
x2=X(2,:);
x3=X(3,:);
% x1=(mapminmax(x1,-1,1));
% x2=(mapminmax(x2,-1,1));
% x3=(mapminmax(x3,-1,1));
save x1 x1;
save x2 x2;
save x3 x3;
save X  X;
save A A;
figure;
subplot(311)
plot(x1);
subplot(312)
plot(x2);
subplot(313)
plot(x3);
% %%=======================================================================%%
N=L;N0=1;
c1_w = LDGT(L,Ls,N,N0,X(1,:),g);
c2_w = LDGT(L,Ls,N,N0,X(2,:),g);
c3_w = LDGT(L,Ls,N,N0,X(3,:),g);
save c1_w c1_w;
save c2_w c2_w;
save c3_w c3_w;
figure;
subplot(131);
contour(abs(c1_w(:,1:N/2)'));
subplot(132);
contour(abs(c2_w(:,1:N/2)'));
subplot(133);
contour(abs(c3_w(:,1:N/2)'));
%------------------------------------------------------------------------%
    
    
    
    
    
    
    
    
    
    
    
    
    
    