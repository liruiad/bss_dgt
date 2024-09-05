function c = DGT(L,N,N0,x,g);
M=L/N0;
c=zeros(M,N);
for m=0:M-1
    RR=zeros(N,1);
    for q=0:N-1
        for p=0:L/N-1
            ii=mod(p*N+q+m*N0,L);
            ix=mod(p*N+q,L);
            RR(q+1)=RR(q+1)+x(ix+1)*g(ii+1);
        end 
    end
c(m+1,:)=fft(RR)';
end
end
