function c = LDGT(L,Ls,N,N0,x,g);
L0=Ls+L;
M=L0/N0;
c=zeros(M,N);
g0=zeros(L0,1);
g0(1:L,1)=g;
x0=zeros(L0,1);
% x0(1:Ls,1)=x;
for m=0:M-1
    RR=zeros(N,1);
    for q=0:N-1
        for p=0:Ls/N-1
            ii=mod(p*N+q-m*N0,L0);
            ix=mod(p*N+q,L0);
            RR(q+1)=RR(q+1)+x(ix+1)*g0(ii+1);
        end 
    end
c(m+1,:)=fft(RR)';
end
end

