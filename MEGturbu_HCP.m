clear all;
path2=[ 'Broad'];
addpath(genpath(path2));

bands = [1 3; 3.5 7; 7.5 13.5; 14 30.5; 31 40];
N=68;

%% create distance matrix from centre of gravity for schaefer1000
load cog_dk68.mat
cog=dk68cog;
for i=1:N
    for j=1:N
        dist(i,j)=sqrt((cog(i,1)-cog(j,1))^2+(cog(i,2)-cog(j,2))^2+(cog(i,3)-cog(j,3))^2);
    end;
end;
[u,v] = find(triu(ones(N),1));    % get edges
% compute distances in same order as ets
for i=1:length(u)
    distance(i)=floor(dist(u(i),v(i)));
end;
[a,idx2]=sort(distance);
dffidx3 = find(diff(a));

%%

% load files
datadir=['Broad/*_MEG_3-Restin*.mat'];
files=dir(datadir);
for ses=1:length(files)
    load(files(ses).name);
    subj{ses}=X;
end

Isubdiag = find(triu(ones(N),-1));
tot=length(Isubdiag);
Isubdiagtot = find(triu(ones(tot),-1));


nn=1;
for f=1:length(bands)
    nn
    flp = bands(f,1);           % lowpass frequency of filter
    fhi = bands(f,2);
    delt = 0.004;            % sampling interval
    k=2;                  % 2nd order butterworth filter
    fnq=1/(2*delt);       % Nyquist frequency
    Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
    [bfilt,afilt]=butter(k,Wn);   % construct the filter
    for ses=1:length(files)
        ses
        HA=double(subj{ses});
        clear signal_filt signal_filt_su signal_filt_syn;
        clear edges edges_su;
        for seed=1:N
            timeserie=filtfilt(bfilt,afilt,detrend(demean(HA(:,seed))));
            signal_filt2=abs(hilbert(zscore(timeserie)));
            signal_filt(seed,:)=zscore(signal_filt2(1000:25:end-1000));
            %             hxs=abs(hilbert(timeserie));
            %             signal_filt2=zscore(movmean(hxs,25));
            %             signal_filt(seed,:)=signal_filt2(1000:25:end-1000);
            Tmaxsu=size(signal_filt,2);
            drawnbin = randi([ceil(Tmaxsu*0.05) ceil(Tmaxsu*0.95)]);
            auxbin = [drawnbin:Tmaxsu 1:drawnbin-1];
            signal_filt_su(seed,:)=signal_filt(seed,auxbin);
            %              signal_filt_su(seed,:)=signal_filt(seed,randperm(size(signal_filt,2)));
        end
        %         signal_filt=zscore(signal_filt,[],2);
        %         signal_filt_su=zscore(signal_filt_su,[],2);
        %         signal_filt_syn=zscore(signal_filt_syn,[],2);
        for t=1:size(signal_filt,2)
            
            %             cphi=cos(signal_filt(:,t));
            %             sphi=sin(signal_filt(:,t));
            %             restacos=repmat(cphi,1,N)-repmat(cphi',N,1);
            %             restasin=repmat(sphi,1,N)-repmat(sphi',N,1);
            %             PhLoMa=sqrt((restacos).^2+(restasin).^2);
            PhLoMa=signal_filt(:,t)*signal_filt(:,t)';
            edges(:,t)=PhLoMa(Isubdiag);
            if ses==length(files)
                stdt_edges(t,:)=mean(PhLoMa);
            end      
            
%             if f==1
%                 save MEG_turbu_HCP_rendering_delta.mat stdt_edges;
%             end
%                         if f==2
%                 save MEG_turbu_HCP_rendering_theta.mat stdt_edges;
%                         end
%                         if f==3
%                 save MEG_turbu_HCP_rendering_alpha.mat stdt_edges;
%                         end
%                         if f==4
%                 save MEG_turbu_HCP_rendering_beta.mat stdt_edges;
%                         end
%                         if f==5
%                 save MEG_turbu_HCP_rendering_gamma.mat stdt_edges;
%             end
%         end
%     end
% end
            
            PhLoMa=signal_filt_su(:,t)*signal_filt_su(:,t)';
            edges_su(:,t)=PhLoMa(Isubdiag);
            if ses==length(files)
                stdt_edges_su(t,:)=mean(PhLoMa);
            end
        end
        if ses==length(files)
            figure;
            maxed=max(max(edges));
            imagesc(edges,[0 maxed]);
            title(['Edges ' num2str(nn)]);
        end        
        if ses==length(files)
            figure;
            imagesc(edges_su,[0 maxed]);
            title(['Edges SURROGATE ' num2str(nn)]);
        end
        
        Edge_meta(nn,ses)=std(edges(:));
        Edge_meta_su(nn,ses)=std(edges_su(:));
        
        gbindist=edges';
        gbindist_su=edges_su';
        ncomm=zeros(tot,tot);
        for k=0:7
            gbindist_k=gbindist(1:end-k,:);
            gbindist_k2=gbindist(1+k:end,:);
            ncomm=corr(gbindist_k,gbindist_k2,'Rows','pairwise');
            ncomm1(k+1,:,:)=ncomm;
            EdgeSpaTimePredictability2(k+1)=nanmean(ncomm(Isubdiagtot));
            
            gbindist_k=gbindist_su(1:end-k,:);
            gbindist_k2=gbindist_su(1+k:end,:);
            ncomm=corr(gbindist_k,gbindist_k2,'Rows','pairwise');
            ncomm1su(k+1,:,:)=ncomm;
            EdgeSpaTimePredictability2_su(k+1)=nanmean(ncomm(Isubdiagtot));
            
        end
        EdgeSpaTimePredictability(nn,ses)=nanmean(EdgeSpaTimePredictability2);
        EdgeSpaTimePredictability_su(nn,ses)=nanmean(EdgeSpaTimePredictability2_su);
        
        if ses==length(files) 
            figure;
            for k=0:7
                for i = 1:length(dffidx3)-1
                    ncomm2{i}=ncomm1(k+1,i,i:end);
                end;
                
                tot2=length(dffidx3);
                for i = 1:tot2
                    tmp=[];
                    for j = 1:tot2-i
                        tmp(j)=ncomm2{j}(i+1);
                    end;
                    avdiagncomm(i)=nanmean(tmp);
                end;
                plot(avdiagncomm,'Linewidth',2); hold on;               
            end;
            title(['Information cascade predictability: ' num2str(nn)]);
            
            figure;
            for k=0:7
                for i = 1:length(dffidx3)-1
                    ncomm2{i}=ncomm1su(k+1,i,i:end);
                end;
                
                tot2=length(dffidx3);
                for i = 1:tot2
                    tmp=[];
                    for j = 1:tot2-i
                        tmp(j)=ncomm2{j}(i+1);
                    end;
                    avdiagncomm(i)=nanmean(tmp);
                end;
                plot(avdiagncomm,'Linewidth',2); hold on;               
            end;
            title(['Information cascade predictability SURROGATES: ' num2str(nn)]);
        end
        
    end
    nn=nn+1;
end

figure
for n=1:length(bands)
    subplot(length(bands),1,n)
    boxplot([Edge_meta(n,:)' Edge_meta_su(n,:)']);
    p(n)=ranksum(Edge_meta(n,:),Edge_meta_su(n,:))
end
sgtitle('Boxplot Edge_meta');


figure
for n=1:length(bands)
    subplot(length(bands),1,n)
    boxplot([EdgeSpaTimePredictability(n,:)' EdgeSpaTimePredictability_su(n,:)']);
    pe(n)=ranksum(EdgeSpaTimePredictability(n,:),EdgeSpaTimePredictability_su(n,:))
end
sgtitle('Boxplot ESP');

save MEG_turbu_HCP.mat stdt_edges_su stdt_edges Edge_meta Edge_meta_su EdgeSpaTimePredictability EdgeSpaTimePredictability_su;

