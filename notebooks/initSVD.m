function [G] = initSVD(Rsum, k)

    [U,D]=eig(Rsum);
    diagDpos=diag(abs(D));
    [Dsort,ind]=sort(diagDpos,1,'descend');%

    trD=sum(diagDpos);
    trDi=0;
    K=1;
    for i=1:length(Dsort)
        trDi = trDi + Dsort(i);
        if trDi/trD > 0.9
            K=i
            break;
        end;
    end;
    k=min(k,K);
    Ginit=[];
    for i=1:k
        xx=U(:,ind(i))*diagDpos(ind(i));
        xp=max(xx,0);
        disp(xp)
        xn=xp-xx;
        if norm(xp)>norm(xn)
            Ginit = [Ginit xp];
        else
            Ginit = [Ginit xn];
        end
    end;
    Gzeros=Ginit<1e-10;
    G=Ginit+1e-10*Gzeros;
end