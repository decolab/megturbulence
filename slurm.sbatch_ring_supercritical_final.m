#!/bin/bash
#SBATCH --qos=vip
#SBATCH --job-name=SuperRing
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-450
#SBATCH --output=Super%A_%a.out
#SBATCH --error=Super%A_%a.err

#Load Matlab 2017a module
ml MATLAB

matlab -nojvm -nodisplay<<-EOF

s=str2num(getenv('SLURM_ARRAY_TASK_ID'))

%% Structure ring

lambda=1;
IntegStepSize=0.1;

NPARCELLS=1000;
Isubdiag = find(triu(ones(NPARCELLS),-1));

rr=zeros(NPARCELLS,NPARCELLS+1);
C=zeros(NPARCELLS,NPARCELLS+1);

for i=1:NPARCELLS
    for j=1:NPARCELLS+1
        d1=abs(j-i);
        d2=NPARCELLS-d1;
        rr(i,j)=min(d1,d2)*IntegStepSize;
    end
end

%% Ordering distances for Edge Matrix
[u,v] = find(triu(ones(NPARCELLS),1));    % get edges
for i=1:length(u)
    distance(i)=floor(2*rr(u(i),v(i)));
end
[adistance,idx2]=sort(distance);
dffidx100 = find(diff(adistance));

for i=1:length(u)
    distance(i)=floor(rr(u(i),v(i)));
end
[adistance,idx2]=sort(distance);
dffidx50 = find(diff(adistance));

%%

for i=1:NPARCELLS
    for j=1:NPARCELLS+1
        C(i,j)=0.5*exp(-lambda*rr(i,j));
    end
end


ni=1;
for i=1:20:NPARCELLS
    nj=1;
    for j=1:20:NPARCELLS
        C50(ni,nj)=mean(C(i,j:j+19));
        nj=nj+1;
    end
    C50(ni,ni)=0.5;
    C50(ni,nj)=C50(ni,1);
    ni=ni+1;
end

ni=1;
for i=1:10:NPARCELLS
    nj=1;
    for j=1:10:NPARCELLS
        C100(ni,nj)=mean(C(i,j:j+9));
        nj=nj+1;
    end
    C100(ni,ni)=0.5;
    C100(ni,nj)=C100(ni,1);
    ni=ni+1;
end

%% Parameters models

beta=2.6;
alpha=atan(beta);
Tmax=100;
dt=0.01;

K=0.05;
omega=ones(NPARCELLS,1);
coupling=K*sqrt(1+beta^2);

NTmax=Tmax/dt+2;

phi=zeros(NPARCELLS,NTmax);
phi100=zeros(100,NTmax);
phi50=zeros(50,NTmax);
Rspatime1000=zeros(NPARCELLS,NTmax);
Rspatime100=zeros(100,NTmax);
Rspatime50=zeros(50,NTmax);
edge_phases100=zeros(length(dffidx100)-1,NTmax);
edge_phases50=zeros(length(dffidx50)-1,NTmax);
phase_lock_matrix_red100=zeros(1,length(dffidx100)-1);
phase_lock_matrix_red50=zeros(1,length(dffidx50)-1);

NSUB=100;

Edge_meta100=zeros(1,NSUB);
Edge_meta50=zeros(1,NSUB);
LocalKoP100=zeros(1,NSUB);
LocalKoP50=zeros(1,NSUB);
LocalKoP1000=zeros(1,NSUB);
KoP1000=zeros(1,NSUB);
StdKoP1000=zeros(1,NSUB);
EdgeSpaTimePredictability100=zeros(NSUB,8);
EdgeSpaTimePredictability50=zeros(NSUB,8);

Sigma_range=0.00001:0.00001:0.0063;
sigma=Sigma_range(s);

D_val=sigma*sqrt(1+beta^2)/K;

for sub=1:NSUB
    sub
    dsig=sqrt(2*sigma*(1+beta^2)*dt);
    phi(:,1)=0.1*ones(NPARCELLS,1);
    for t=0:dt:600
        phiaux=vertcat(phi(:,1),phi(1,1));
        phiaux2=complex(cos(phiaux),-sin(phiaux));
        phimat=repmat(phiaux2',NPARCELLS,1);
        LK=trapz(C.*phimat,2)*IntegStepSize;
        R=abs(LK);
        Theta=wrapTo2Pi(angle(LK));
        phi(:,1)=wrapTo2Pi(phi(:,1)+dt*(omega-coupling*R.*sin(phi(:,1)-Theta+alpha))+dsig*randn(NPARCELLS,1));
    end
    nn=1;
    for t=0:dt:Tmax
        phiaux=vertcat(phi(:,nn),phi(1,nn));
        phiaux2=complex(cos(phiaux),-sin(phiaux));
        phimat=repmat(phiaux2',NPARCELLS,1);
        LK=trapz(C.*phimat,2)*IntegStepSize;
        R=abs(LK);
        Theta=wrapTo2Pi(angle(LK));
        phi(:,nn+1)=wrapTo2Pi(phi(:,nn)+dt*(omega-coupling*R.*sin(phi(:,nn)-Theta+alpha))+dsig*randn(NPARCELLS,1));
        nn=nn+1;
    end
    
    %%%
    for i=1:NPARCELLS
        if i<=NPARCELLS/2
            fracG(i)=1/(1+(2*pi*(i-1)/NPARCELLS/IntegStepSize)^2);
        else
            fracG(i)=1/(1+(2*pi*(NPARCELLS-1-i)/NPARCELLS/IntegStepSize)^2);
        end
    end
    fracG=fracG';
    
    nn=1;
    for t=0:dt:Tmax
        phiaux=complex(cos(phi(:,nn)),sin(phi(:,nn)));
        LK=ifft(IntegStepSize*fracG.*fft(phiaux))/IntegStepSize;
        R=abs(LK);
        Rspatime1000(:,nn)=abs(LK);
        Theta=wrapTo2Pi(angle(LK));
        phi(:,nn+1)=wrapTo2Pi(phi(:,nn)+dt*(omega-coupling*R.*sin(phi(:,nn)-Theta+alpha))+dsig*randn(NPARCELLS,1));
        nn=nn+1;
    end
    %%%
    
    nn=1;
    for n=1:10:NPARCELLS   
        phi100(nn,:)=wrapTo2Pi(atan2(nanmean(sin(phi(n:n+9,:))),nanmean(cos(phi(n:n+9,:)))));
        nn=nn+1;
    end

    nn=1;
    for n=1:20:NPARCELLS
        phi50(nn,:)=wrapTo2Pi(atan2(nanmean(sin(phi(n:n+19,:))),nanmean(cos(phi(n:n+19,:)))));
        nn=nn+1;
    end
    
    %%%
    Rglobal=abs(sum(complex(cos(phi),sin(phi)),1))/NPARCELLS;
    KoP1000(sub)=mean(Rglobal);
    StdKoP1000(sub)=std(Rglobal);
    
    for nn=1:NTmax
        phiaux=vertcat(phi(:,nn),phi(1,nn));
        phiaux2=complex(cos(phiaux),-sin(phiaux));
        phimat=repmat(phiaux2',NPARCELLS,1);
        LK=trapz(C.*phimat,2)*IntegStepSize;
        Rspatime1000(:,nn)=abs(LK);
        
        phiaux=vertcat(phi100(:,nn),phi100(1,nn));
        phiaux2=complex(cos(phiaux),-sin(phiaux));
        phimat=repmat(phiaux2',100,1);
        LK=trapz(C100.*phimat,2)*IntegStepSize;
        Rspatime100(:,nn)=abs(LK);
        
        phiaux=vertcat(phi50(:,nn),phi50(1,nn));
        phiaux2=complex(cos(phiaux),-sin(phiaux));
        phimat=repmat(phiaux2',50,1);
        LK=trapz(C50.*phimat,2)*IntegStepSize;
        Rspatime50(:,nn)=abs(LK);
        
        %%% Edge-centric 100
        
        cphi=cos(phi(:,nn));
        sphi=sin(phi(:,nn));
        restacos=repmat(cphi,1,NPARCELLS)-repmat(cphi',NPARCELLS,1);
        restasin=repmat(sphi,1,NPARCELLS)-repmat(sphi',NPARCELLS,1);
        PhLoMa=sqrt((restacos).^2+(restasin).^2);
        phase_lock_matrix=PhLoMa(Isubdiag);
        for i = 1:length(dffidx100)-1
            phase_lock_matrix_red100(i)=nanmean(phase_lock_matrix(dffidx100(i)+1:dffidx100(i+1)));
        end
        edge_phases100(:,nn)=phase_lock_matrix_red100;
        
        %%% Edge-centric 50
        
        cphi=cos(phi(:,nn));
        sphi=sin(phi(:,nn));
        restacos=repmat(cphi,1,NPARCELLS)-repmat(cphi',NPARCELLS,1);
        restasin=repmat(sphi,1,NPARCELLS)-repmat(sphi',NPARCELLS,1);
        PhLoMa=sqrt((restacos).^2+(restasin).^2);
        phase_lock_matrix=PhLoMa(Isubdiag);
        for i = 1:length(dffidx50)-1
            phase_lock_matrix_red50(i)=nanmean(phase_lock_matrix(dffidx50(i)+1:dffidx50(i+1)));
        end
        edge_phases50(:,nn)=phase_lock_matrix_red50;
    end
    LocalKoP1000(sub)=std(Rspatime1000(:));
    LocalKoP100(sub)=std(Rspatime100(:));
    LocalKoP50(sub)=std(Rspatime50(:));
    
    %%% Edge-centric 100
    
    Edge_meta100(sub)=std(edge_phases100(:));
    
    tot=length(dffidx100)-1;
    gbindist=edge_phases100';
    ncomm=zeros(tot,tot);
    for k=0:7
        for i = 1:tot
            for j = 1:tot
                ccaux=corrcoef(gbindist(1:end-k,i),gbindist(1+k:end,j),'Rows','pairwise');
                ncomm(i,j)=ccaux(2);
            end
        end
        EdgeSpaTimePredictability100(sub,k+1)=nanmean(nanmean(ncomm));
    end
    
    %%% Edge-centric 50
    
    Edge_meta50(sub)=std(edge_phases50(:));
    
    tot=length(dffidx50)-1;
    gbindist=edge_phases50';
    ncomm=zeros(tot,tot);
    for k=0:7
        for i = 1:tot
            for j = 1:tot
                ccaux=corrcoef(gbindist(1:end-k,i),gbindist(1+k:end,j),'Rows','pairwise');
                ncomm(i,j)=ccaux(2);
            end
        end
        EdgeSpaTimePredictability50(sub,k+1)=nanmean(nanmean(ncomm));
    end    
end

save(sprintf('Wsuper_%03d.mat',s),'EdgeSpaTimePredictability100','EdgeSpaTimePredictability50','Edge_meta100','Edge_meta50','LocalKoP100','LocalKoP50','LocalKoP1000','KoP1000','StdKoP1000','D_val');
EOF


