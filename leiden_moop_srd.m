%% Leiden ranking based demonstration of PCA, standard and sparse Non-negative matrix factorization and SRD
% written by 
% Uses nnmfv1_4 - The non-negative matrix factorization toolbox for biological data mining


close all
clear all
clc

%% Adding nnmfv1_4 toolbox

folder = fileparts(which('leiden_moop_srd.m')); 
addpath(genpath(folder));

%% Loading Data
[numread,textread] = xlsread('CWTS Leiden Ranking 2020_v1_2.xlsx','Sorted'); %Leiden 2020 
unames = textread(2:end,1);
unames = strrep(unames,"University", "");
unames = strrep(unames,"of", "");
unames = strrep(unames,"  ", " ");
varnames=textread(1,3:end);
nat = textread(2:end,2);
X=numread;
%Entropy
EntropyR = [];
HyperVolume = [];
%% Scale to [0;1]
X(isnan(X))=0;
N=size(X,1);
n=size(X,2)


X=(X-repmat(min(X),N,1))./repmat((max(X)-min(X)),N,1);

%% Calculate the ranks; let's use tiderank
% required by the nonparametric tests signrank and ranksum, and for the computation of Spearman's rank correlation
R=[];
for k=1:n
    R=[R tiedrank(-X(:,k))];
end    

%% Problem Formulation
%Distribution of the whole set.
xx=-X;
p_rank=zeros(N,1);
index=1:N;
l=0;
while min(p_rank)==0
    l=l+1;
    mem=paretosort(xx);
    imem=find(mem==1);
    p_rank(index(imem))=l;
    xx(imem,:)=[];
    index(imem)=[];
end
[Dr,r]= hist(p_rank,unique(p_rank));
index=find(p_rank==l); 
EntropyR=[EntropyR sum(Dr/N.*log10(Dr/N))/log10(1/N)]; 
HyperVolume =[HyperVolume hypervolume(X(index, :), max(X), 10000)];


% Draw RE figure - first subplot
entrFig = figure(1);
ax = gca;
ax.Visible = 'off';
subplot(2, 3, 1);
ax = bar(r,Dr, 'LineWidth',0.01);
ax.EdgeColor = 'none';
t = title("Distribution of the set;" + " RE(D) = " + num2str(EntropyR(1), 4))
t.FontSize = 20;
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
%Calculating the objective pair with the highest entropy.
v = 1:n;
C = nchoosek(v,2);
parvar = [];
RE = [];
%RE = sum([X(:,C(1,1)) X(:, C(1, 2))]/N.*log10([X(:, C(1,1)) X(:, C(1, 2))]/N))/log10(1/N);
for i = 1:length(C)
%RE = [RE sum(((R(:,C(i,1)))/N).*log10(R(:, C(i,1)))/N)/sum(log10(R(:, C(i, 2)))/N)];
xx=[X(:, C(i,1)), X(:, C(i,2))];
p_rank=zeros(N,1);
index=1:N;
l=0;
while min(p_rank)==0
    l=l+1;
    mem=paretosort(-xx);
    imem=find(mem==1);
    p_rank(index(imem))=l;
    xx(imem,:)=[];
    index(imem)=[];
end
[Dr,r]= hist(p_rank,unique(p_rank));
RE=[RE sum(Dr/N.*log10(Dr/N))/log10(1/N)];
end

[RE, relabels] = sort(RE);
EntropyR = [EntropyR RE(end)];
parvar = C(relabels(end), 1:2);
figure(2)
handle = scatter(X(:,parvar(1)), X(:,parvar(2)), 10, 'filled');
axis([0 1.2 0 1.2])
xlabel("Total Number of Open Access publications, " + varnames{parvar(1)}, 'Interpreter', 'none');
ylabel("Number of Green Access publcations, " + varnames{parvar(2)}, 'Interpreter', 'none');
%'linewidth', 15,'Color', [205 100 100]/255
xx=[X(:, parvar(1)), X(:, parvar(2))];
p_rank=zeros(N,1);
index=1:N;
l=0;
while min(p_rank)==0
    l=l+1;
    mem=paretosort(-xx);
    imem=find(mem==1);
    p_rank(index(imem))=l;
    xx(imem,:)=[];
    index(imem)=[];
end
nl=4
colors={'r','k'};
for l=1:nl 
    index=find(p_rank==l); 
    
    if(l == 1)
        HyperVolume =[HyperVolume hypervolume([X(index, parvar(1)), X(index, parvar(2))], max([X(index, parvar(1)), X(index, parvar(2))]), 10000)]; 
    end
    
    hold on
    for k=1:length(index) 
        if(l == 1)
        text(X(index(k),parvar(1)),X(index(k),parvar(2)),unames{index(k)},'Color',colors{mod(l,2)+1}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize',20);
        scatter(X(index(k),parvar(1)), X(index(k),parvar(2)),30,'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        %set(handle, 'MarkerIndices', markers(1, index(k)) , 'MarkerSize',15, 'Color', colors{l})
        elseif(l == 2)
        text(X(index(k),parvar(1)),X(index(k),parvar(2)),unames{index(k)},'Color',colors{mod(l,2)+1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize',20);
        scatter(X(index(k),parvar(1)), X(index(k),parvar(2)),30, 'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor', [1 0 0]);
        elseif(l == 3)
        text(X(index(k),parvar(1)),X(index(k),parvar(2)),unames{index(k)},'Color',colors{mod(l,2)+1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize',20)
        scatter(X(index(k),parvar(1)), X(index(k),parvar(2)),30, 'MarkerFaceColor',[0 0 0], 'MarkerEdgeColor', [0 0 0]);
        else
        text(X(index(k),parvar(1)),X(index(k),parvar(2)),unames{index(k)},'Color',colors{mod(l,2)+1}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize',20) ;
        scatter(X(index(k),parvar(1)), X(index(k),parvar(2)),30, 'MarkerFaceColor', [1 0 0] , 'MarkerEdgeColor', [1 0 0]);
        end
    end    
end
%Drawing lines & Paretos for the two objectives
% A = [0.014 0.209];
% B = [0.008 0.179]
% xlim = get(gca,'XLim');
% m = (B(2)- B(1))/(A(2)-A(1));
% n = B(2)*m - A(2);
% y1 = m*xlim(1) + n;
% y2 = m*xlim(2) + n;
% hold on
% line([xlim(1) xlim(2)],[y1 y2])

ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
figure(1)
subplot(2, 3, 6);
[Dr,r]= hist(p_rank,unique(p_rank));
ax = bar(r,Dr, 'LineWidth',0.01);
ax.EdgeColor = 'none';
t = title("Highest Entropy pair;" + " RE(D) = " + num2str(EntropyR(2),4))
t.FontSize = 20;
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;

%genDistance = [genDistance, generational_distance([X(:, parvar(1)), X(:, parvar(2))])];



%% TOPSIS

A =X'; 
A = normalize(A, 'norm');
best = max(A);
worst = min(A);
Dmax = sqrt(sum((A - repmat(best, size(A, 1), 1)).^2)); %best ideal
Dmin = sqrt(sum((A - repmat(worst, size(A, 1), 1)).^2)); %worst ideal
[tps, tpslabels] = sort(Dmax./(Dmin + Dmax)); 
[tpsw, tpslabelsw] = sort(Dmin./(Dmin + Dmax));


%
figure(10);
scatter(1:N, tps, 16,[0 0 0],'filled');
hold on
scatter(1:10, tps(1, 1:10), 30,[1 0 0],'filled');
scatter((N-10):N, tps(1, (end-10):end), 30,[0 0 1],'filled');
xlabel('University Rank');
ylabel('L_2 Distance to best condition');
axis([1 N 0.3 0.8])
for i=1:10
    text(150, 0.3 + i*0.03, unames{tpslabels(1, i)}, 'FontSize',30);
    quiver(i, tps(i),  (150-i),(0.3 + i*0.03)- tps(i), 'MaxHeadSize', 0, 'Color', 'k')
end

xticks([1 100 200 300 400 500 600 700 800 900 1000 1100 1200])
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;

for i=1:10
    
   unames{tpslabelsw(1, i)};
  
end
%Comment no tied values
tpsSpearman = 1 - (6* sum((tiedrank(tps) -(1:N)).^2)./N*(N^2-1)); 

%% SRD - Objectives 
gmax= max(numread,[],2);
figure(20)
[nrow,ncol]=size(R);
%max srd
    if rem(nrow,2)==1
        k=(nrow-1)/2;
        m=2*k*(k+1);
    else
        k=nrow/2;
        m=2*k^2;
    end
%the best "virtiual" method  / ideal objective / ideal ranking  
nrk=tiedrank(-gmax); 
names=varnames;
% Calculate the SRD
srd=sum(abs(R-repmat(nrk,1,n)),1)/m*100;
[srdi,si]=sort(srd);
nsrdi=names(si); 
srdi;
% Simulation of the probabolity distribution 
nSim=1e4;
S=[];
io=[1:nrow];
for i=1:nSim
       S(i)=sum(abs(io-randperm(nrow)))/m*100;
end
%yyaxis left
[prob,srdc]= histcounts(S,unique(S),'Normalization','cdf');
plot([srdc],[0 prob] * 100)
xlabel('SRD')
ylabel('P(SRD) [%]')
hold on 
% Calculate the probabilities of the orderings (based on cum prob)
psrdi = interp1([0 srdc],[ 0 0 prob], srdi);
psrdi = psrdi +1;
yyaxis right
axv=axis;
ylabel('SRD')
%axis([0 2*max(srdi) 0 1.5*max(psrdi)])
for i=1:n 
    
   % text(srdi(i),srdi(i),nsrdi{i},'fontsize',11)
line([srdi(i) srdi(i) ],[0 srdi(i)])
hold on 
end
%out=srdi;

%% Removing the worst objectives - G7-9 C10 O13 (column 34, 41-43)

remX = [X(:, 1:33) X(:,35:40) X(:, 44:46)]
redVarLabels = {varnames{1, 1:33} varnames{1, 35:40} varnames{1, 44:46}}

%% PCA 
%Apply PCA
[COEFF, SCORE, latent, tsquared, explained] = pca(zscore(numread));

%Draw PCA components
figure(30)
bar(explained, 'EdgeColor',[0 .5 .5],'LineWidth',0.01);
xlabel('Component Number');
ylabel('Variance [%]');

%Draw PCA biplot
figure(31) 
handles = biplot([COEFF(:, 1), COEFF(:, 2)],'Scores', [SCORE(:, 1) SCORE(:, 2)],'VarLabels', varnames);
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize', 14, 'fontWeight', 'bold');
axis([-0.08 0.30 -0.25 0.28]);
xlabel('PC 1 (51.96%)');
ylabel('PC 2 (13.91%)');
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;


figure(1202)% only explains between 50-60%
scatter(SCORE(:, 1), SCORE(:, 2), 10, 'filled');
xlabel('PC 1 (51.96%)');
ylabel('PC 2 (13.91%)');
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;


xx2=-[SCORE(:, 1) SCORE(:, 2)];
p_rank=zeros(N,1);
index=1:N;
l=0;
while min(p_rank)==0
    l=l+1;
    mem=paretosort(xx2);
    imem=find(mem==1);
    p_rank(index(imem))=l;
    xx2(imem,:)=[];
    index(imem)=[];
end
hold on 
nl=2
colors={'r', 'k'}

scPar = (SCORE-repmat(min(SCORE),N,1))./repmat((max(SCORE)-min(SCORE)),N,1);

for l=1:nl
    hold on 
    index=find(p_rank==l); 
    [linV linLbl] = sort(SCORE(index, 1));
    if(l == 1)
    HyperVolume = [HyperVolume hypervolume(scPar(index, 1:2), max(scPar(:, 1:2)), 10000)];
    end
    plot(SCORE(index(linLbl), 1), SCORE(index(linLbl), 2), 'k--');
    for k=1:length(index)
        if(l== 1)
        an = text(SCORE(index(k), 1), SCORE(index(k), 2), unames{index(k)},'Color',colors{mod(l,2)+1}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize',20);
        set(an,'Rotation',45);
        scatter(SCORE(index(k), 1), SCORE(index(k), 2), 30,'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        else
        an = text(SCORE(index(k), 1), SCORE(index(k), 2), unames{index(k)} ,'Color', colors{mod(l,2)+1}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize',20); %num2str(index(k))
        set(an,'Rotation',45);  
        scatter(SCORE(index(k), 1), SCORE(index(k), 2), 30,'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
        end
    end    
end
annotation('textarrow',[0.3/3.5 0.4245/3.5], [2/2.5 1.8668/2.5],'Linewidth',2)

%Draw to entropy figure
figure(1) %SCORE(:, 2)
subplot(2, 3, 4);
[Dr,r]= hist(p_rank, unique(p_rank));
EntropyR=[EntropyR sum(Dr/N.*log10(Dr/N))/log10(1/N)]%%
ax = bar(r,Dr, 'LineWidth',0.01);
ax.EdgeColor = 'none';
t = title("PCA;" + " RE(D) = " + num2str(EntropyR(3)), 4);
t.FontSize = 20;
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;

%PCA heatmap - Supplementary
figure(200)
heatmap({'PC1 (51.96%)', 'PC2 (13.91%)', 'PC3 (8.49%)', 'PC4 (5.36%)', 'PC5 (4.45%)', 'PC6 (2.58%)', 'PC7 (2.23%)'}, varnames, round(COEFF(:, 1:7), 4), 'Colormap', pink);
ax = gca;
ax.FontSize = 16;

 %% NNMF - Scaled Data
rng (4) % Setting random for Reproductibility
opt = statset('maxiter',5000, 'display','final');

[Y2,H2] = nnmf(X,2,'opt',opt,'alg','als');

%% Scaled NNMF Figures
%biplot
figure(50);
ych = find(contains(nat,'China'));
yjap = find(contains(nat,'Japan'));
yus = find(contains(nat,'United States'));
yko = find(contains(nat,'Korea'));
yge = find(contains(nat,'Germany'));
yir = find(contains(nat,'Iran'));
yta = find(contains(nat,'Taiwan'));
yrem = 1:N
yrem = yrem'
yrem(ych) = 0
yrem(yjap) = 0
yrem(yus) = 0
yrem(yko) = 0
yrem(yge) = 0
yrem(yir) = 0
yrem(yta) = 0

yrem(find(yrem == 0)) = []

plot(Y2(yrem,1),Y2(yrem,2), '.', 'MarkerSize',25, 'MarkerFaceColor', [0 0.4470 0.7410]);
hold on
plot(Y2(ych,1),Y2(ych,2), 'd', 'MarkerSize',10, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor',[0 0 0]);
plot(Y2(yjap,1),Y2(yjap,2), 'o','MarkerSize',10, 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerEdgeColor',[0 0 0]);
plot(Y2(yus,1),Y2(yus,2),'p','MarkerSize',15, 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'MarkerEdgeColor',[0 0 0]);
plot(Y2(yko,1),Y2(yko,2),'*','MarkerSize',15, 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
plot(Y2(yge,1),Y2(yge,2),'^','MarkerSize',15, 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerEdgeColor',[0 0 0]);
plot(Y2(yir,1),Y2(yir,2),'s','MarkerSize',15, 'MarkerFaceColor', [0 0.5 1], 'MarkerEdgeColor',[0 0 0]);
plot(Y2(yta,1),Y2(yta,2),'h','MarkerSize',15, 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerEdgeColor',[0 0 0]);


xlabel("w_1");
ylabel("w_2");


xx=-Y2;
p_rank=zeros(N,1);
index=1:N;
l=0;
while min(p_rank)==0
    l=l+1;
    mem=paretosort(xx);
    imem=find(mem==1);
    p_rank(index(imem))=l;
    xx(imem,:)=[];
    index(imem)=[];
end
hold on 
nl=1
colors={'r','k'}
%ures = []
for l=1:nl
    hold on 
    index=find(p_rank==l);
    HyperVolume = [HyperVolume hypervolume(Y2(index, :), max(Y2), 10000)]
    [linV linLbl] = sort(Y2(index, 1));
    plot(Y2(index(linLbl), 1), Y2(index(linLbl), 2), 'k--');
    for k=1:length(index)
        an = text(Y2(index(k),1),Y2(index(k),2),unames{index(k)},'Color',colors{2},  'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize',20) ;% num2str(index(k)) unames{index(k)}  
        %scatter(Y2(index(k), 1), Y2(index(k), 2), 60,'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        set(an,'Rotation',45);
        
    end    
end

annotation('textarrow',[0.3/3.5 0.4245/3.5], [2/2.5 1.8668/2.5],'Linewidth',2)

%quiver(0.3, 2, 0.4245 -0.3, 1.8668-2, 'k', 'MaxHeadSize', 10, 'lineWidth', 5)

[ax icons] = legend({'Unspecified','China', 'Japan', 'United States', 'South Korea', 'Germany', 'Iran', 'Taiwan', '1st Pareto front'}, 'FontSize',18);
icons = findobj(icons,'Type','patch');
icons = findobj(icons,'Marker','none','-xor');
set(icons,'MarkerSize', 10);

D = norm(X - Y2*H2,'fro')/sqrt(N*n) %residual
axis([0 3.5 0 2.5]);
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;

figure(51)
biplot(H2','scores',Y2,'varlabels',varnames);
axis([0 0.6 0 0.6]);
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',6);

figure(1)
subplot(2, 3, 2);
[Dr,r]= hist(p_rank,unique(p_rank));
EntropyR=[EntropyR sum(Dr/N.*log10(Dr/N))/log10(1/N)];
ax = bar(r,Dr, 'LineWidth',0.01);
ax.EdgeColor = 'none';
t = title("NNMF with Scaled Input;" + " RE(D) = " + num2str(EntropyR(4), 4));
t.FontSize = 20;
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;

%% NNMF - Removed Objectives

opt = statset('maxiter', 5000, 'display','final');
[remY, remH] = nnmf(remX,2,'opt',opt,'alg','als');

figure(52)
scatter(remY(:,1),remY(:,2),10, [0 0 0], 'filled');
hold on
scatter(remY(ych,1),remY(ych,2),30,'filled', 'd');
scatter(remY(yjap,1),remY(yjap,2),30,'o', 'LineWidth', 2);
scatter(remY(yus,1),remY(yus,2),60, '+', 'LineWidth', 2);
scatter(remY(yko,1),remY(yko,2),60,'*');
scatter(remY(yge,1),remY(yge,2),30, 'x', 'LineWidth',6);
scatter(remY(yir,1),remY(yir,2),30, 'filled','s');
scatter(remY(yta,1),remY(yta,2),60,[1 0 0], 'h','filled', 'LineWidth', 6);
xlabel('W_1');
ylabel('W_2');

xx=-remY;
p_rank=zeros(N,1);
index=1:N;
l=0;
while min(p_rank)==0
    l=l+1;
    mem=paretosort(xx);
    imem=find(mem==1);
    p_rank(index(imem))=l;
    xx(imem,:)=[];
    index(imem)=[];
end
hold on 
nl=1
colors={'r','k'}
%ures = []
for l=1:nl
    hold on 
    index=find(p_rank==l) ;
    [linV linLbl] = sort(remY(index, 1));
    plot(remY(index(linLbl), 1), remY(index(linLbl), 2), 'k--');
    for k=1:length(index)
        an = text(remY(index(k),1),remY(index(k),2),unames{index(k)},'Color',colors{2},  'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize',20) ;% num2str(index(k)) unames{index(k)}  
        scatter(remY(index(k), 1), remY(index(k), 2), 60,'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        set(an,'Rotation',45);
        
    end    
end

[ax icons] = legend({'Remaining', 'China', 'Japan', 'United States', 'South Korea', 'Germany', 'Iran', 'Taiwan', '1st Pareto front'}, 'FontSize',18);
icons = findobj(icons,'Type','patch');
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',20);
axis([0 3 0 7]);
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;

%% NNMF - Ranked
%opt = statset('maxiter',200,'display','final');
%[w,h] = nnmf(R,2,'rep',25,'opt',opt,'alg','mult');
opt = statset('maxiter',5000,'display','final');
[Y,H,residualr] = nnmf(R,2,'alg','als');
figure(100)
scatter(Y(:,1),Y(:,2),10, [0 0 0], 'filled');

figure(101)
biplot(H','scores',Y,'varlabels',varnames);
axis([0 0.35 0 0.35])
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize', 14, 'fontWeight', 'bold');
xlabel("Component 1");
ylabel("Component 2");
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;

figure(1)
subplot(2, 3, 3);
[Dr,r]= hist(p_rank,unique(p_rank));
EntropyR=[EntropyR sum(Dr/N.*log10(Dr/N))/log10(1/N)];
ax = bar(r,Dr, 'LineWidth',0.01);
ax.EdgeColor = 'none';
t = title("NNMF with Ranked Input;" + " RE(D) = " + num2str(EntropyR(5), 4));
t.FontSize = 20;

for l=1:nl
    hold on 
    index=find(p_rank==l) ;
    HyperVolume = [HyperVolume hypervolume(Y(index, :), min(Y), 10000)]
end

ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;

%% Sparse with removed objectives
 %ANNM - W | YNNM - H

option.eta =  0.4; % 0.4
option.beta = 0.4; % 4
option.dis=true;
option.residual=1e-14;
option.tof=1e-14;
[sW,sH,numIter,tElapsed,finalResidual] = sparsenmfnnls(remX, 2, option);

sprValue = [sparsity(sW), sparsity(sH)];
figure(600)
plot(sW(yrem,2),sW(yrem,1), '.', 'MarkerSize',25, 'MarkerFaceColor', [0 0.4470 0.7410]);
hold on
plot(sW(ych,2),sW(ych,1), 'd', 'MarkerSize',10, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor',[0 0 0]);
plot(sW(yjap,2),sW(yjap,1), 'o','MarkerSize',10, 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerEdgeColor',[0 0 0]);
plot(sW(yus,2),sW(yus,1),'p','MarkerSize',15, 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'MarkerEdgeColor',[0 0 0]);
plot(sW(yko,2),sW(yko,1),'*','MarkerSize',15, 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
plot(sW(yge,2),sW(yge,1),'^','MarkerSize',15, 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerEdgeColor',[0 0 0]);
plot(sW(yir,2),sW(yir,1),'s','MarkerSize',15, 'MarkerFaceColor', [0 0.5 1], 'MarkerEdgeColor',[0 0 0]);
plot(sW(yta,2),sW(yta,1),'h','MarkerSize',15, 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerEdgeColor',[0 0 0]);
xlabel('First Weight Component');
ylabel('Second Weight Component');

xx=-sW;
p_rank=zeros(N,1);
index=1:N;
l=0;
while min(p_rank)==0
    l=l+1;
    mem=paretosort(xx);
    imem=find(mem==1);
    p_rank(index(imem))=l;
    xx(imem,:)=[];
    index(imem)=[];
end
hold on 
nl=1
colors={'r','k'}
%ures = []
for l=1:nl
    hold on 
    index=find(p_rank==l) ;
    HyperVolume = [HyperVolume hypervolume(sW(index, :), max(sW), 10000)]
    [linV linLbl] = sort(sW(index, 2));
    plot(sW(index(linLbl), 2), sW(index(linLbl), 1), 'k--');
    for k=1:length(index)
        an = text(sW(index(k),2),sW(index(k), 1),unames{index(k)},'Color',colors{2},  'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize',20) ;% num2str(index(k)) unames{index(k)}  
        scatter(sW(index(k), 2), sW(index(k), 1), 60,'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        set(an,'Rotation',45);
        
    end    
end

[ax icons] = legend({'Unspecified', 'China', 'Japan', 'United States', 'South Korea', 'Germany', 'Iran', 'Taiwan', '1st Pareto front'}, 'FontSize',18);
icons = findobj(icons,'Type','patch');
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',20);
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;

% Sparse biplot
figure(601)
biplot(sH','scores', sW,'varlabels',redVarLabels);
axis([0 3.3 0 3.3])
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize', 10, 'fontWeight', 'bold');
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;



figure(1) %SCORE(:, 2)
subplot(2, 3, 5)
[Dr,r]= hist(p_rank,unique(p_rank));
EntropyR=[EntropyR sum(Dr/N.*log10(Dr/N))/log10(1/N)]
%genDistance = [genDistance, generational_distance(sW2)];
ax = bar(r,Dr, 'LineWidth',0.01);
ax.EdgeColor = 'none';
t = title("SNNMF with Scaled Input;" + " RE(D) = " + num2str(EntropyR(6)), 4)
t.FontSize = 20;
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;

figure(1)
ax = axes(entrFig);
ax.Visible = 'off';
ax.XLabel.Visible = 'on';
ax.YLabel.Visible = 'on';
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
ylabel('D(r) - number of samples with rank r');
xlabel('r - rank');

sH = [sH(2, :); sH(1, :)]
figure(602)
%heatmap(redVarLabels,{'H1', 'H2'}, sH);
heatmap(redVarLabels,{'H1', 'H2'},  round(sH, 2));
ax = gca;
ax.FontSize = 16;


figure(700)

scatter(EntropyR, HyperVolume, 40, [0 0 0], 'filled')
ylabel('Hypervolume');
xlabel('Relative Entropy');

annotation('textarrow',[0.1 0.1], [0.3 0.1] ,'String',' ','FontSize',13,'Linewidth',2)
%annotation('textarrow',[0.1 0.3], [0.05 0.05] ,'String',' ','FontSize',13,'Linewidth',2)
methods = {' Set', ' Highest Entropy pair', ' PCA', ' NNMF - scaled', ' NNMF - ranked', ' SNNMF'}

for i=1:6
    text(EntropyR(i), HyperVolume(i), methods{i}, 'FontSize', 20);
    %quiver(i, tps(i),  (150-i),(0.3 + i*0.03)- tps(i), 'MaxHeadSize', 0, 'Color', 'k')
end

ax = gca
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
