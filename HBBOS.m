function [u0,yy,xa]=HBBOS (f,a,b,d,MAXGEN,N)
%混合蛙跳和生物地理学融合的混合优化算法
m=5 ;%种群分组数
n=fix(N/m);%每组青蛙包含的个数
N=m*n;
Ne=n; %组内迭代数
dn=d+1;
X=rand(N,d+1);
aa=repmat(a,N,1);
bb=repmat(b,N,1);
X(:,1:d)=aa+(bb-aa).*X(:,1:d);
X(:,d+1)=feval(f, X(:,1:d));
X = PopSort(X);

yy=zeros(1,MAXGEN);
lambda=(1:n)/n;
lambdaMin = min(lambda);
lambdaMax = max(lambda);
lambdaScale = (lambda - lambdaMin) / (lambdaMax - lambdaMin);

for ii=1:MAXGEN   
    Ped=1-(1-0.96)*ii/MAXGEN;    
    for i4=1:m %分组         
        local = X(i4:m:end,:);
        Island=local(:,1:d);        
        for k=1:Ne %每组青蛙迭代次数
            if k==1
                alfa1=rand;
                rnum=ceil(n*rand);while rnum==k,rnum=ceil(n*rand);end
                rnum1=ceil(n*rand);while rnum1==k ||rnum1==rnum,rnum1=ceil(n*rand);end
                Island(k,:) = local(k,1:d)+alfa1*(X(1,1:d)-local(k,1:d)+local(rnum1,1:d)-local(rnum,1:d));
            else
                for j=1:d
                    SelectIndex=ceil((k-1)*rand);Num=k;
                    if rand < lambdaScale(k)
                        % Pick a habitat from which to obtain a feature                        
                        alfa=rand;                        
                        if rand>Ped%i>M/2
                            cnum=ceil(rand*d);
                            Island(k,j) =alfa*local(SelectIndex,j)+(1-alfa)*local(SelectIndex,cnum);                            
                        else
                            Island(k,j) =local(SelectIndex,j)+2*(alfa-0.5)*(local(SelectIndex,j)-local(k,j));
                        end                       
                    else                        
                        if ii>MAXGEN/2
                            Island(k,j) =local(k,j)+rand*(X(1,j)-local(k,j)+local(SelectIndex,j)-local(Num,j));                            
                        else
                            Num=ceil(n*rand);while Num==SelectIndex,Num=ceil(n*rand);end
                            Island(k,j) =local(k,j)+2*(rand-0.5)*(X(1,j)-local(k,j)+local(SelectIndex,j)-local(Num,j));                            
                        end                           
                    end
                end
            end
        end%结束Ne
        Island=ControlBoundD(Island,aa(1:n,:),bb(1:n,:));        
        fit=feval(f,Island);
        for k=1:Ne
            value=fit(k);
            if local(k,dn)>value
                local(k,dn)=value;
                local(k,1:d)=Island(k,:);
            end
        end
        X(i4:m:end,:) =local;
    end   %结束m    
    X = PopSort(X);
    yy(ii)=X(1,dn);
end     %结束MAXGEN
xa=X(1,1:d);u0=X(1,dn);
