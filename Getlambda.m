function [Idx,lambda_reduced,lh1,lh2,lh3,rat11,rat22,rat33]=Getlambda(lambda1,locs,d1,d2,d3,n);
[Idx,C]=kmeans(lambda1',n);
%Idx=[3 3 1 2 1 1 1];
C=C';
[row col]=size(lambda1);
Ck=n;
posv=zeros(Ck,1);
lambda2=zeros(2,Ck);
lh1=zeros(2,col); i1=1;rat11=zeros(2,col);
lh2=zeros(2,col); i2=1;rat22=zeros(2,col);
lh3=zeros(2,col); i3=1;rat33=zeros(2,col);
for ix=1:Ck
    v=[];
    Lhv0=0;
    Lhv1=0;
    lambda2(:,ix)=0;
    for iy=1:col 
          if(Idx(iy)==ix)
             v=[v iy];
             a_j0=lambda1(1,iy);
             a_j1=lambda1(2,iy);
             Lhv00=(1/(sqrt(2*pi)*1)*exp(-1/(2*1)*abs(a_j0*d1(locs(iy))-d2(locs(iy)))^2/(1+a_j0^2)));
             Lhv11=(1/(sqrt(2*pi)*1)*exp(-1/(2*1)*abs(a_j1*d1(locs(iy))-d3(locs(iy)))^2/(1+a_j1^2)));
             if(Lhv00>Lhv0)
                lambda2(1,ix) = lambda1(1,iy);
                Lhv0=Lhv00;
             end
             if(Lhv11>Lhv1)
                lambda2(2,ix) = lambda1(2,iy);
                Lhv1=Lhv11;
             end
             
             if(ix==1)
                 lh1(1,i1)= Lhv00;
                 lh1(2,i1)= Lhv11;
                 rat11(1,i1)=a_j0;
                 rat11(2,i1)=a_j1;
                 i1=i1+1;
             end
             if(ix==2)
                 lh2(1,i2)= Lhv00;
                 lh2(2,i2)= Lhv11;
                 rat22(1,i2)=a_j0;
                 rat22(2,i2)=a_j1;
                 i2=i2+1;
             end
             if(ix==3)
                 lh3(1,i3)= Lhv00;
                 lh3(2,i3)= Lhv11;
                 rat33(1,i3)=a_j0;
                 rat33(2,i3)=a_j1;
                 i3=i3+1;
             end  
          end
    end
    posv(ix,1) = min(v);
end
lambda = zeros(2,Ck);
sposv=sort(posv);
posv1=zeros(Ck,1);
for ix=1:Ck;
    posv1(ix)=find(sposv== posv(ix));
end
for ix=1:Ck;
    iy=find(posv1==ix);
    lambda(:,ix)=lambda2(:,iy);
end
lambda_reduced=lambda2;
end