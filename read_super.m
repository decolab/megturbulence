clear all;

for s=1:450
    fileName = sprintf('Wsuper_%03d.mat',s);
    if exist(fileName, 'file') == 2
        load(fileName);
        Edge_meta100_all(s,:)=Edge_meta100;
        EdgeSpaTimePredictability100_all(s,:,:)=EdgeSpaTimePredictability100;
        KoP1000_all(s,:)=KoP1000;
        StdKoP1000_all(s,:)=StdKoP1000;
        LocalKoP1000_all(s,:)=LocalKoP1000;
        LocalKoP100_all(s,:)=LocalKoP100;
        D_range(s)=D_val;
    else
        Edge_meta100_all(s,:)=NaN*ones(1,100);
        EdgeSpaTimePredictability100_all(s,:,:)=NaN*ones(100,8);
        KoP1000_all(s,:)=NaN*ones(1,100);
        StdKoP1000_all(s,:)=NaN*ones(1,100);
        LocalKoP1000_all(s,:)=NaN*ones(1,100);
        LocalKoP100_all(s,:)=NaN*ones(1,100);
        D_range(s)=NaN;
    end
end

EdgeSpaTimePredictability100_all1=nanmean(EdgeSpaTimePredictability100_all,3);

figure(1)
shadedErrorBar(D_range,nanmean(LocalKoP1000_all,2),nanstd(LocalKoP1000_all,[],2),'-k',0.7)
hold on;
shadedErrorBar(D_range,nanmean(KoP1000_all,2),nanstd(KoP1000_all,[],2),'-b',0.7)
shadedErrorBar(D_range,nanmean(Edge_meta100_all,2),nanstd(Edge_meta100_all,[],2),'-r',0.7)
shadedErrorBar(D_range,nanmean(StdKoP1000_all,2),nanstd(StdKoP1000_all,[],2),'-m',0.7)
shadedErrorBar(D_range,40*nanmean(EdgeSpaTimePredictability100_all1,2),nanstd(EdgeSpaTimePredictability100_all1,[],2),'-g',0.7)

figure(2)
shadedErrorBar(D_range,nanmean(LocalKoP100_all,2),nanstd(LocalKoP100_all,[],2),'-k',0.7)
hold on;
shadedErrorBar(D_range,nanmean(Edge_meta100_all,2),nanstd(Edge_meta100_all,[],2),'-r',0.7)
shadedErrorBar(D_range,nanmean(LocalKoP1000_all,2),nanstd(LocalKoP1000_all,[],2),'-g',0.7)


%% Tests . limites 0.025 -- 0.18
i=1;
refD=324; %% 324
for dd=D_range
   [pl1000(i) h_LKop1000ref(i)]=ranksum(LocalKoP1000_all(refD,:),LocalKoP1000_all(i,:),'Tail','left','Alpha',0.05);
   [pl100(i) h_LKop100ref(i)]=ranksum(LocalKoP100_all(refD,:),LocalKoP100_all(i,:),'Tail','left','Alpha',0.05);
   [pe100(i) h_Edge100ref(i)]=ranksum(Edge_meta100_all(refD,:),Edge_meta100_all(i,:),'Tail','left','Alpha',0.05);
   [p h_ESP100ref(i)]=ranksum(EdgeSpaTimePredictability100_all1(refD,:),EdgeSpaTimePredictability100_all1(i,:),'Tail','left','Alpha',0.05);
   [p h_ESP100ref2(i)]=ranksum(EdgeSpaTimePredictability100_all1(1,:),EdgeSpaTimePredictability100_all1(i,:),'Tail','right','Alpha',0.05);
   i=i+1;
end

figure(4)
subplot(2,1,1);
plot(D_range,h_LKop1000ref,'-m');
hold on;
plot(D_range,h_Edge100ref*2,'-r');
subplot(2,1,2);
plot(D_range,h_LKop100ref,'-m');
hold on;
plot(D_range,h_Edge100ref*2,'-r');
% plot(D_range,h_ESP100ref*2,'-g');

figure(5)
esp100=squeeze(nanmean(EdgeSpaTimePredictability100_all,2));
imagesc(D_range,1:8,esp100');

figure(6)
shadedErrorBar(D_range,nanmean(EdgeSpaTimePredictability100_all1,2),nanstd(EdgeSpaTimePredictability100_all1,[],2),'-g',0.7)
ranksum(EdgeSpaTimePredictability100_all1(1,:),EdgeSpaTimePredictability100_all1(49,:))
ranksum(EdgeSpaTimePredictability100_all1(49,:),EdgeSpaTimePredictability100_all1(322,:))
ranksum(EdgeSpaTimePredictability100_all1(322,:),EdgeSpaTimePredictability100_all1(end,:))

save results_super.mat D_range Edge_meta100_all EdgeSpaTimePredictability100_all KoP1000_all StdKoP1000_all LocalKoP1000_all LocalKoP100_all;


%%% Figures paper
%% BOXPLOTS;

figure(1)
subplot(4,1,1);
boxplot([KoP1000_all(1,:)' KoP1000_all(94,:)' KoP1000_all(300,:)' KoP1000_all(324,:)']);
ylim([0 1.2]);
subplot(4,1,2);
boxplot([StdKoP1000_all(1,:)' StdKoP1000_all(94,:)' StdKoP1000_all(300,:)' StdKoP1000_all(324,:)']);
ylim([0 0.03]);
subplot(4,1,3);
boxplot([LocalKoP1000_all(1,:)' LocalKoP1000_all(94,:)' LocalKoP1000_all(300,:)' LocalKoP1000_all(324,:)']);
ylim([0 0.3]);
subplot(4,1,4);
boxplot([LocalKoP100_all(1,:)' LocalKoP100_all(94,:)' LocalKoP100_all(300,:)' LocalKoP100_all(324,:)']);
ylim([0 0.3]);

ranksum(KoP1000_all(300,:),KoP1000_all(325,:))
ranksum(StdKoP1000_all(300,:),StdKoP1000_all(325,:))

figure(2)
boxplot([Edge_meta100_all(1,:)' Edge_meta100_all(94,:)' Edge_meta100_all(300,:)' Edge_meta100_all(324,:)']);
ylim([0 0.75]);
figure(3)
boxplot([Edge_meta100_all(94,:)' Edge_meta100_all(300,:)' Edge_meta100_all(324,:)']);

figure(4)
boxplot([EdgeSpaTimePredictability100_all1(1,:)' EdgeSpaTimePredictability100_all1(94,:)' EdgeSpaTimePredictability100_all1(300,:)' EdgeSpaTimePredictability100_all1(end,:)']);

%% Corr turbu esp edgemeta

n=1;
for i=90:323
    edgy(n)=nanmean(LocalKoP1000_all(i,:));
    espy(n)=nanmean(EdgeSpaTimePredictability100_all1(i,:));
    n=n+1;
end

figure(5)
scatter(edgy,espy);
    