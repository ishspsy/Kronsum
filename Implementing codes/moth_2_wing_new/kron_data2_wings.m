
addpath('umich_file')



Total_Data=cell(2,7);

Total_Data{1,1}=csvread('TphaseensembleJ.csv');
Total_Data{1,2}=csvread('TphaseensembleK.csv');
Total_Data{1,3}=csvread('TphaseensembleL.csv');
Total_Data{1,4}=csvread('TphaseensembleM.csv');
Total_Data{1,5}=csvread('TphaseensembleN.csv');
Total_Data{1,6}=csvread('TphaseensembleP.csv');
Total_Data{1,7}=csvread('TphaseensembleQ.csv');

Total_Data{2,1}=csvread('TspiketrigensembleJ.csv');
Total_Data{2,2}=csvread('TspiketrigensembleK.csv');
Total_Data{2,3}=csvread('TspiketrigensembleL.csv');
Total_Data{2,4}=csvread('TspiketrigensembleM.csv');
Total_Data{2,5}=csvread('TspiketrigensembleN.csv');
Total_Data{2,6}=csvread('TspiketrigensembleP.csv');
Total_Data{2,7}=csvread('TspiketrigensembleQ.csv');


Y_Data{1,1}=csvread('wingstroke_data_J.csv');
Y_Data{1,2}=csvread('wingstroke_data_K.csv');
Y_Data{1,3}=csvread('wingstroke_data_L.csv');
Y_Data{1,4}=csvread('wingstroke_data_M.csv');
Y_Data{1,5}=csvread('wingstroke_data_N.csv');
Y_Data{1,6}=csvread('wingstroke_data_P.csv');
Y_Data{1,7}=csvread('wingstroke_data_Q.csv');

Y_Data{2,1}=csvread('wingstroke_data_J.csv');
Y_Data{2,2}=csvread('wingstroke_data_K.csv');
Y_Data{2,3}=csvread('wingstroke_data_L.csv');
Y_Data{2,4}=csvread('wingstroke_data_M.csv');
Y_Data{2,5}=csvread('wingstroke_data_N.csv');
Y_Data{2,6}=csvread('wingstroke_data_P.csv');
Y_Data{2,7}=csvread('wingstroke_data_Q.csv');


for ii=1:14;
rdata=Total_Data{ii};
X=rdata; Y=Y_Data{ii}; Y=[Y(:,end),Y(:,end-1)];
[n,p]=size(X);
ind=find(sum(isnan(X),2));
X=X(setdiff(1:n,ind),:); Y=Y(setdiff(1:n,ind),:); [n,p]=size(X); Y_dif=Y(:,1)-Y(:,2); indd=setdiff(1:n,find(isnan(Y_dif)));
X=X(indd,:); Y=Y(indd,:); [n,p]=size(X); Y_dif=Y(:,1)-Y(:,2); Y_dif2=[(1:n)',Y_dif]; Y_dif2=sortrows(Y_dif2,-2); ind=Y_dif2(:,1);
X=X(ind,:); Y=Y(ind,:);
X=X-mean(mean(X));
Total_Datao{ii}=X;
X=X-repmat(mean(X,2),[1,p]);
mmm=mean(mean(X.^2)); mmm=sqrt(mmm);
X=sqrt(2)*X./mmm;
%X=X-repmat(mean(X,2),[1,p]);
Total_Data{ii}=X; Y_Data{ii}=Y;
end;
%save('Total_Dat3.mat', 'Total_Datao','Total_Data', 'Y_Data')

load('Total_Dat3.mat') 
load('fin_moth_gp_structure_for_Rsquared2.mat')



%%
load('Total_Dat3.mat')
load('fin_moth_gp_structure_for_Rsquared2.mat')
%ii=2;  ii=6; ii=4; ii=14;

X=Total_Data{ii}; ttt=corr(X');
clf;
figure 
imagesc(ttt); 
colormap(b2r(-1, 1)), colorbar;
tt=title('Moth J: spatial-correlation');
xlabel('Wingstroke'); ylabel('Wingstroke')
set(tt, 'FontSize', 18);
print -depsc Moth_14_spat_corr


X=Total_Data{ii}; ttt=corr(X);
clf;
figure 
imagesc(ttt); 
colormap(b2r(-1, 1)), colorbar;
tt=title('Moth J: temporal-correlation');
xlabel('Time'); ylabel('Time')
set(tt, 'FontSize', 18);
print -depsc Moth_6_temp_corr



%%

title1=cell(1,16);
title1{1}='Moth 1'; title1{2}='Moth 1';
title1{3}='Moth 2'; title1{4}='Moth 2';
title1{5}='Moth 3'; title1{6}='Moth 3';
title1{7}='Moth 4'; title1{8}='Moth 4';
title1{9}='Moth 5'; title1{10}='Moth 5';
title1{11}='Moth 6'; title1{12}='Moth 6';
title1{13}='Moth 7'; title1{14}='Moth 7';
title1{15}='General model'; title1{16}='General model';

%graphical structure

tauA_set=[1,1.3,1.5,1.7,1.9]; lamA_set=[0.0001 0.0003 0.0005 0.001 0.003 0.005 0.01 0.02];
for ii=[1,2,5,6];
Gp_B_set=cell(1,length(tauA_set));  Gp_IB_set=cell(1,length(tauA_set));
Gp_A_set=cell(1,length(tauA_set));  Gp_IA_set=cell(1,length(tauA_set));
Gp_A_set1=cell(1,length(lamA_set));Gp_IA_set1=Gp_A_set1;Gp_B_set1=Gp_A_set1;Gp_IB_set1=Gp_A_set1;    
for jj=1:length(tauA_set);
    parfor kk=1:length(lamA_set)
XX=Total_Data{ii}; lamA=lamA_set(kk);
tauA=tauA_set(jj); tauB=2-tauA; [n,p]=size(XX); ini_X=XX'*XX/n-tauB*eye(p); eta=0.5*norm(ini_X); 
[G_A, G_IA]=estimator_i2_lassoo(XX,tauA,tauB,eta,lamA,1);
ini_X2=XX*XX'/p-tauA*eye(n); eta=0.5*norm(ini_X2); 
[G_B, G_IB]=estimator_i2_lassoo(XX',tauB,tauA,eta,lamA,1);
Gp_A_set1{kk}=G_A; Gp_IA_set1{kk}=G_IA; Gp_B_set1{kk}=G_B; Gp_IB_set1{kk}=G_IB;       
    end
Gp_A_set{jj}=Gp_A_set1; Gp_IA_set{jj}=Gp_IA_set1; Gp_B_set{jj}=Gp_B_set1;Gp_IB_set{jj}=Gp_IB_set1;
end
GGp_A_set{ii}=Gp_A_set; GGp_IA_set{ii}=Gp_IA_set;
GGp_B_set{ii}=Gp_B_set; GGp_IB_set{ii}=Gp_IB_set;
end

%save('fin_moth_gp_structure_for_Rsquared2_1256.mat','GGp_A_set','GGp_IA_set','GGp_B_set','GGp_IB_set')
%load('fin_moth_gp_structure_for_Rsquared2.mat')
%load('fin_moth_gp_structure_for_Rsquared2_1256.mat')

for ii=1:14
    for jj=1:5
        GGp_A_set{ii}{jj}=(GGp_A_set{ii}{jj}+GGp_A_set{ii}{jj}')*0.5;
        GGp_B_set{ii}{jj}=(GGp_B_set{ii}{jj}+GGp_B_set{ii}{jj}')*0.5;
        if min(eig(GGp_A_set{ii}{jj}))<0
            GGp_A_set{ii}{jj}=GGp_A_set{ii}{jj}-(min(eig(GGp_A_set{ii}{jj}))-0.1)*eye(500);
        end
        if min(eig(GGp_B_set{ii}{jj}))<0
            GGp_B_set{ii}{jj}=GGp_B_set{ii}{jj}-(min(eig(GGp_B_set{ii}{jj}))-0.1)*eye(size(GGp_B_set{ii}{jj},2));
        end       
        GGp_A_set{ii}{jj}=GGp_A_set{ii}{jj}*500*tauA_set(jj)/trace(GGp_A_set{ii}{jj});
        GGp_IA_set{ii}{jj}=inv(GGp_A_set{ii}{jj});
        GGp_B_set{ii}{jj}=GGp_B_set{ii}{jj}*size(GGp_B_set{ii}{jj},2)*(2-tauA_set(jj))/trace(GGp_B_set{ii}{jj});
        GGp_IB_set{ii}{jj}=inv(GGp_B_set{ii}{jj});
    end
end



%%%%%

load('Moth_5_graph_A.mat'); ii=5   % n=298; p=500;
load('Moth_1_graph_A.mat'); ii=1   % n=535; p=500;

%aaa=importdata('fin_moth_gp_structure_for_Rsquared2.mat')
aaa=importdata('fin_moth_gp_structure_for_Rsquared2_1256.mat')
tauA_set=[1,1.3,1.5,1.7,1.9]; lamA_set=[0.0001 0.0003 0.0005 0.001 0.003 0.005 0.01 0.02];

%G_IA=aaa.GGp_IA_set{6}{1};  G_IA=G_IA-diag(diag(G_IA)); mmax=max(max((G_IA)));mmin=min(min((G_IA)));
G_IA=aaa.GGp_IA_set{2}{3}{6}; G_IA=full(G_IA); G_IA=G_IA-diag(diag(G_IA)); %mmax=max(max((G_IA)));mmin=min(min((G_IA)));
G_IA=aaa.GGp_IA_set{6}{3}{6}; G_IA=full(G_IA); G_IA=G_IA-diag(diag(G_IA)); %mmax=max(max((G_IA)));mmin=min(min((G_IA)));
%G_IA=G_IA.*(abs(G_IA)>0.05);
%quantile(abs(offidiag1),0.4)



clf;
figure 
imagesc(G_IA); 
colormap(b2r(mmin,mmax)), colorbar;
tt=title('Moth J: $$\hat{\Theta}$$','Interpreter','Latex');
%tt=title(sprintf( ' %s: correlation matrix' ,  title1{i}),'Interpreter','Latex');
xlabel('Time'); ylabel('Time')
set(tt, 'FontSize', 18);
print -depsc Moth_1_IA


%
G_IA=G_IA-diag(diag(G_IA)); B=G_IA; n=size(G_IA,1);
for i=1:(n-1);
    %B(i,i+1)=0.5;
    %B(i+1,i)=0.5;
    B(i,i+1)=4.7;
    B(i+1,i)=4.7;
end

iindx=[]; iindy=[];
for i=1:n;
for j=1:i;
%if abs(B(i,j))>0.25; %i=1
%if abs(B(i,j))>0.1; %i=5
if abs(B(i,j))>0.3; %i=5
iindx=[iindx,i]; iindy=[iindy,j];
else;
iindx=iindx; iindy=iindy;
end;
end;
end;
real_ind=union(iindx,iindy);

B=B(real_ind,real_ind);
n=size(B,1);
iindx=[]; iindy=[];
for i=1:n;
for j=1:i;
if  abs(B(i,j))>0.3;     %abs(B(i,j))>2.7;
iindx=[iindx,i]; iindy=[iindy,j];
else;
iindx=iindx; iindy=iindy;
end;
end;
end;

clf
G=graph(iindx,iindy);
h = plot(G,'Layout','force');
%labelnode(h, 1:size(real_ind,2),real_ind)
h.XData=3*h.XData;
h.YData=3*h.YData;
h.NodeLabel=real_ind;
h.MarkerSize=2
set(gca,'fontsize',2)
set(gca,'Visible','off')
%print -depsc Moth_1_IB
print -depsc Moth_3_IA_force2
%




load('Mooth_1_graph_B.mat'); 
load('Mooth_5_graph_B.mat'); 

B=G_IB-diag(diag(G_IB)); B=B.*(abs(B)>0.15); i=1; y_d=Y_Data{i};
B=G_IB-diag(diag(G_IB)); B=B.*(abs(B)>0.01); i=1; y_d=Y_Data{i};  %second
%B=G_IB-diag(diag(G_IB)); B=B.*(abs(B)>0.25); i=1; y_d=Y_Data{i};
B=G_IB-diag(diag(G_IB)); B=B.*(abs(B)>0.05); i=5; y_d=Y_Data{i};


G_IB=aaa.GGp_IB_set{2}{3}{7}; G_IB=full(G_IB); G_IB=G_IB-diag(diag(G_IB)); %mmax=max(max((G_IA)));mmin=min(min((G_IA)));
G_IB=aaa.GGp_IB_set{6}{3}{6}; G_IB=full(G_IB); G_IB=G_IB-diag(diag(G_IB)); %mmax=max(max((G_IA)));mmin=min(min((G_IA)));


B=G_IB-diag(diag(G_IB)); n=size(B,1);
n=size(B,1);  ng=round(n/3); ind1=1:ng; ind2=ng+1:2*ng; ind3=2*ng+1:n; AB=abs(B);
[sum(sum(AB(ind1,ind1))),sum(sum(AB(ind2,ind2))),sum(sum(AB(ind3,ind3)))]
[sum(sum(AB(ind1,ind2))),sum(sum(AB(ind1,ind3))),sum(sum(AB(ind2,ind3)))]


B=G_IB-diag(diag(G_IB)); n=size(B,1);
iindx=[]; iindy=[];
for i=1:n;
for j=1:i;
%if abs(B(i,j))>0.25; %i=1
if abs(B(i,j))>0.15; %i=5
iindx=[iindx,i]; iindy=[iindy,j];
else;
iindx=iindx; iindy=iindy;
end;
end;
end;
real_ind=union(iindx,iindy);



B=B(real_ind,real_ind);
n=size(B,1);
iindx=[]; iindy=[];
for i=1:n;
for j=1:i;
if abs(B(i,j))>0.15;
iindx=[iindx,i]; iindy=[iindy,j];
else;
iindx=iindx; iindy=iindy;
end;
end;
end;

clf
G=graph(iindx,iindy);
h = plot(G,'Layout','force');
%labelnode(h, 1:size(real_ind,2),real_ind)
h.NodeLabel=real_ind;
h.MarkerSize=3
set(gca,'fontsize',10)
set(gca,'Visible','off')
%print -depsc Moth_1_IB
print -depsc Moth_5_IB


bins = conncomp(G);
bins2=[real_ind',bins']; 
bin_fin=sortrows(bins2,2); bin_fin=[bin_fin(:,2), bin_fin(:,1)];
bin_num=max(bin_fin(:,1)); bin_siz=[];
for i=1:bin_num
bin_siz=[bin_siz, sum(bin_fin(:,1)==i)];
end
bin_sizz=[(1:bin_num)',bin_siz'];  ord_bin=sortrows(bin_sizz,-2); 

%moth 1
i=1; i=2; i=4
ind=find(bin_fin(:,1)==ord_bin(i,1)); clu_gp=bin_fin(ind,2)
i=3;i=7
i=5; i=9

%moth 5
ind=find(bin_fin(:,1)==ord_bin(i,1)); clu_gp=bin_fin(ind,2); clu_gp'



%moth1
I1=[39,44, 53,107, 112, 117,140, 148,159]; I1=[39 44 112 159]; I1=[107 117 140 148]
I2=[227, 268,283,300,332,389];   I2=[161 227 283]; I2=[268 332 389]
I3=[411, 432, 483, 487,525]; I3=[432 483 487]
I4=[246 39 221 112 353 307 44 159];
I5=[75 58 127 161 300 227 283];

%new moth1
I1=[58 75 127 161 36 488 134 302];
I2=[137 44 159  166  167 8];
I3=[413 509 321 428 327];

I1=[238 404 445 399 359 310  331 377 193  399 313 391 341 435];
I2=[303 171 48 354 270 120 398 51];
I3=[161 127 58 75  488 36 134 302 515];

I1=[415 44 443 54 191 158 82 186 113 114 112 135 122 53 192];
I2=[28 389 414 459 2 283 392 461 140 181 336 193 398 11 130 155 213 354 509 102 74 325 527 234 31 71 189];
I3=[255 388  275 343 273];

X=Total_Datao{2};  
X1=X(I1,:); X2=X(I2,:); X3=X(I3,:); X1=X1*10^4;  X2=X2*10^4;  X3=X3*10^4; 
X4=X(I4,:); X5=X(I5,:);
YD=Y_Data{1}(:,1)-Y_Data{1}(:,2); [mean(YD(I1)), mean(YD(I2)), mean(YD(I3))], [std(YD(I1)), std(YD(I2)), std(YD(I3))]

%moth 5
I1=[ 12    18    44    69    83   102 ];
I3=[ 208   218   227   237   248   253   263   287];
I2=[104   129   144];
I4=[26  73 80 91 92 185 216];

I1=[10 29 44 14  50 16 15 51  12 31 60 158 175 22 140 161 52 129 27 42];
I2=[297 281 295 39 98 107 4 171  128 162];
I3=[9 139 26 18  233 58 78];

X=Total_Data{6};
X1=X(I1,:); X2=X(I2,:); X3=X(I3,:); X4=X(I4,:); X1=X1*10^4;  X2=X2*10^4;  X3=X3*10^4; 
YD=Y_Data{5}(:,1)-Y_Data{5}(:,2); YD=10*YD; [mean(YD(I1)), mean(YD(I2)), mean(YD(I3))], [std(YD(I1)), std(YD(I2)), std(YD(I3))]




clf
subplot(1,3,1)
h1=plot(X1','b--', 'Linewidth', 0.2); hold on;
xlim([1 500]); %ylim([-4 4])
set(gca,'FontSize', 15)
subplot(1,3,2)
h2=plot(X2','k--','Linewidth', 0.2); hold on
xlim([1 500]); %ylim([-4 4])
xlabel('Time') %,'FontSize',   23) % x-axis label
ylabel('Torque')
tt=title('Moth 3: Wingstrokes'); tt.FontSize=13;
set(gca,'FontSize', 15)
subplot(1,3,3)
h3=plot(X3','r--','Linewidth', 0.2);  %hold on;
xlim([1 500]); %ylim([-4 4])
set(gca,'FontSize', 15)
%h4=plot(X4','c-.','Linewidth', 1); hold on;

%print -depsc moth_2_wing_new2
%print -depsc moth_1_wing
print -depsc moth_6_wing_new2









clf;
figure 
imagesc(G_IB-diag(diag(G_IB))); 
colormap(b2r(-0.1, 0.1)), colorbar;
tt=title('Moth J: $$\hat{\Theta}$$','Interpreter','Latex');
%tt=title(sprintf( ' %s: correlation matrix' ,  title1{i}),'Interpreter','Latex');
set(tt, 'FontSize', 18);
print -depsc Moth_1_IB


%%%


kk=3;
%spl_mat=zeros(500,length(-20:3:520)-kk-1);
spl_mat=zeros(500,length(-200:25:700)-kk-1);
for x=1:500;
for j=1:(length(-200:25:700)-kk-1)
[vy,vx] = bspline_basis(j,kk,-200:25:700,x); %jntx
spl_mat(x,j)=vy;
end;
end;

spl_mat=spl_mat(:,6:27);
%csvwrite('splmat2.csv',spl_mat);

%spl_mat=spl_mat(:,3:end-2);   %csvwrite('splmat.csv',spl_mat);  

spl_mat=csvread('splmat2.csv');
%spl_mat=[spl_mat, 0.01*(1:500)', 0.0001*((1:500).^2)',0.000001*((1:500).^3)'];

[YYY,R] = qr(spl_mat,0); 
spl_mat=YYY;



[U,S,V] = svd(spl_mat,'econ');
CC=V*S*S*V';
C=sqrtm(CC); tauC=trace(C)/size(C,1);
projj=spl_mat*inv(spl_mat'*spl_mat)*spl_mat';



tauA=1; X_var=cell(1,14);
for ii=1:14;
aaa11=Y_Data{ii}; [n,p]=size(aaa11);
X=Total_Data{ii};  var_X=mean(X.^2,1)-mean(X,1).^2; 
y=aaa11;  y2=y(:,1)-y(:,2);  [y2]=simp_lin_reg(Total_Datao{ii}, y2);
y2=2*sqrt(size(y2,1))*y2/norm(y2);
XX=X*spl_mat;    var_XX=mean(XX.^2,1)-mean(XX,1).^2; 

lamr=0:3:80; 
parfor jr=1:length(lamr)
bb=ridge(y2,XX,lamr(jr)); b=inv(XX'*XX+ lamr(jr)*eye(size(XX,2)))*XX'*y2; 
b=spl_mat*b;  %ii=1
beta_ridge_d{jr}=b;  %[norm(bbeta_ridge_d{ii}),norm(y2-X*bbeta_ridge_d{ii})]
Rs2_ridge(jr)=corr(X*b,y2)^2;
end
bbeta_ridge_d{ii}=beta_ridge_d;
Rsq2_ridge{ii}=Rs2_ridge;
laml=0.01:0.2:4; laml=[0, 0.000001, 0.000005, 0.00001, 0.00005,0.0001,0.0005,0.001,0.005,laml,5,6,7,8,9,10,11, 12,13,14,15, 16,20,30,50,80,100];
parfor jl=1:length(laml)
[B,FitInfo] = lasso(XX,y2,'Lambda', laml(jl)*sqrt(log(size(XX,2))/size(XX,1))); b2=spl_mat*B;
beta_lasso_d{jl}=b2;
Rs2_lasso(jl)=corr(X*b2,y2)^2;
end
bbeta_lasso_d{ii}=beta_lasso_d;
Rsq2_lasso{ii}=Rs2_lasso;
end



tauA_set=0.1:0.1:1.9;

for ii=1:12
Gp_B_set=cell(1,length(tauA_set));  Gp_IB_set=cell(1,length(tauA_set));
Gp_A_set=cell(1,length(tauA_set));  Gp_IA_set=cell(1,length(tauA_set));
    
parfor jj=1:length(tauA_set);
XX=Total_Data{ii}; %XX=XX-repmat(mean(XX),[size(XX,1),1]);
tauA=tauA_set(jj); tauB=2-tauA; [n,p]=size(XX); ini_X=XX'*XX/n-tauB*eye(p); eta=0.5*norm(ini_X); %eta=eta_set(ii);
[G_A, G_IA]=estimator_i2_lassoo(XX,tauA,tauB,eta,0.01,1);
ini_X2=XX*XX'/p-tauA*eye(n); eta=0.5*norm(ini_X2);  %eta=eta_set(ii);
[G_B, G_IB]=estimator_i2_lassoo(XX',tauB,tauA,eta,0.01,1);
if min(eig(G_IA))<0
    G_IA=G_IA-(min(eig(G_IA))-0.3)*eye(size(G_IA,1));
end
if min(eig(G_IB))<0
    G_IB=G_IB-(min(eig(G_IB))-0.3)*eye(size(G_IB,1));
end

G_A=inv(G_IA); G_A=G_A*tauA*p/trace(G_A); G_IA=inv(G_A);
G_B=inv(G_IB); G_B=G_B*tauB*n/trace(G_B); G_IB=inv(G_B);

Gp_A_set{jj}=G_A; Gp_IA_set{jj}=G_IA; Gp_B_set{jj}=G_B; Gp_IB_set{jj}=G_IB;
end
GGp_A_set{ii}=Gp_A_set; GGp_IA_set{ii}=Gp_IA_set;
GGp_B_set{ii}=Gp_B_set; GGp_IB_set{ii}=Gp_IB_set;
end




const_den_set_d=[1 1 1 1 1 1 1 1 1 1 1 1 1 1]; const_num_set_d=0.5*[2 2 2 2 2 2 2 2 2 2 2 2 2 2];
lam_set=[0.0001, 0.0005, 0.001,0.01,0.1];

for ii=1:12;
tauA_set=[1,1.3,1.5,1.7,1.9]; 
trace_A_set=zeros(1,length(tauA_set)); min_eig_A_set=zeros(1,length(tauA_set));
max_eig_A_set=zeros(1,length(tauA_set)); Rsq2_EIV_set=zeros(1,length(tauA_set));Rsq2_EIV_set2=zeros(1,length(tauA_set));
EIV_coef=cell(1,length(tauA_set));
parfor jj=1:length(tauA_set);
const_den_d=const_den_set_d(ii);
const_num_d=const_num_set_d(ii);
XX=Total_Data{ii}; 
tauA=tauA_set(jj); tauB=2-tauA; [n,p]=size(XX); ini_X=XX'*XX/n-tauB*eye(p); eta=0.5*norm(ini_X); %eta=eta_set(ii);
G_A=GGp_A_set{ii}{jj};G_IA=GGp_IA_set{ii}{jj}; G_A=(G_A+G_A')/2;G_IA=(G_IA+G_IA')/2;
G_B=GGp_B_set{ii}{jj};G_IB=GGp_IB_set{ii}{jj}; G_B=(G_B+G_B')/2;G_IB=(G_IB+G_IB')/2;
aaa11=Y_Data{ii}; [n,p]=size(aaa11);
X=Total_Data{ii}; y=aaa11;  
y2=y(:,1)-y(:,2);  [y2]=simp_lin_reg(Total_Datao{ii}, y2);
y2=2*sqrt(size(y2,1))*y2/norm(y2);
XX=X*spl_mat; 
tauACB=(norm(XX,'fro')^2)/(size(XX,1)*size(XX,2)); tau_A_tilde=tauACB-tauB*tauC;
[n,p]=size(XX); Gamma=XX'*XX/n-tauB*C; eta=0.1*norm(Gamma);
%[est_beta,proc]=composite_regression(tau_A_tilde,tauB,C,XX,y2,eta,0.01,1,randn(p,1)); bet=spl_mat*est_beta;
[est_beta,proc]=composite_regression(tau_A_tilde,tauB,C,XX,y2,eta,0.2,1,randn(p,1)); bet=spl_mat*est_beta;
est_beta=bet;  EIV_coef{jj}=est_beta; y=y2;
[rsq1, rsq2]=R_squared_calculation(X,y,G_IA,G_IB,est_beta,n,tauB,const_den_d,const_num_d);
Rsq2_EIV_set(jj)=rsq1; Rsq2_EIV_set2(jj)=rsq2;


end
RRsq2_EIV_set{ii}=Rsq2_EIV_set;RRsq2_EIV_set2{ii}=Rsq2_EIV_set2; EEIV_coef{ii}=EIV_coef;
end


%%
%%
denomi=sum((y-mean(y)).^2)*n*(1/min(eig(G_IA)))*norm(est_beta)^2;
numera_x=1+p*log(n)/(4*n); numera_x=1;
numera_x2=1+length(find(est_beta))*log(n)/(4*n);
denomi2=numera_x*sum((y-mean(y)).^2)*n*(1/min(eig(G_IA)))*norm(est_beta)^2;     
numera=max(0,abs(y'*X*est_beta)-abs(sqrt(numera_x2)*norm(y)*sqrt(n)*sqrt(1/min(eig(G_IB)))*norm(est_beta)));    
numera2=max(0,abs(y'*X*est_beta)-const*abs(norm(y)*sqrt(n)*sqrt(1/min(eig(G_IB)))*norm(est_beta)));
den_true=sum((y-mean(y)).^2)*est_beta'*Xp'*Xp*est_beta;
EIV_RS22=[EIV_RS22,numera2^2/denomi2];

%%


ind_las=[19,14,23,7,22,7];
for ii=1:12
clf
%h1=plot(1:500, 10*bbeta_ridge_d{ii}{7},'k-', 'Linewidth', 1); hold on;
h1=plot(1:500, 10*bbeta_ridge_d{ii}{7},'k-', 'Linewidth', 1); hold on;

%h2=plot(1:500, bbeta_lasso_d{ii}{ind_las(round(ii/2))},'b--', 'Linewidth', 1); hold on;
h2=plot(1:500, 10*bbeta_lasso_d{ii}{13},'b--', 'Linewidth', 1); hold on;
h3=plot(1:500, 1*EEIV_coef{ii}{5}, 'r-.', 'Linewidth', 1); hold on
xlim([1 500]); %ylim([-4 4])
xla=xlabel('Time') %,'FontSize',   23) % x-axis label
yla=ylabel('Value')
xla.FontSize=12; yla.FontSize=12;
%tt=title('Moth P: Coefficient','FontSize',15);
tt=title(sprintf('%s: Coefficient',title1{ii}),'Interpreter','Latex');
tt.FontSize=15;
hlegend=legend('Ridge', 'Lasso' , 'EIV', 'Location','north');  hlegend.FontSize=15;
legend('boxoff')
print([sprintf('newmoths003_coef%d_3',ii )],'-depsc')
%print([sprintf('newmoths002_coef%d_3',ii )],'-depsc')
%print -depsc moth_12_coef
%print([sprintf('moths_coef%d',ii )],'-depsc')
end




load('Mooth_1_graph_B.mat'); 
load('Mooth_5_graph_B.mat'); 

%B=G_IB-diag(diag(G_IB)); B=B.*(abs(B)>0.25); i=1; y_d=Y_Data{i};
B=G_IB-diag(diag(G_IB)); B=B.*(abs(B)>0.15); i=5; y_d=Y_Data{i};

n=size(B,1);  ng=round(n/3); ind1=1:ng; ind2=ng+1:2*ng; ind3=2*ng+1:n; AB=abs(B);
[sum(sum(AB(ind1,ind1))),sum(sum(AB(ind2,ind2))),sum(sum(AB(ind3,ind3)))]
[sum(sum(AB(ind1,ind2))),sum(sum(AB(ind1,ind3))),sum(sum(AB(ind2,ind3)))]


B=G_IB-diag(diag(G_IB));
iindx=[]; iindy=[];
for i=1:n;
for j=1:i;
%if abs(B(i,j))>0.25; %i=1
if abs(B(i,j))>0.25; %i=5
iindx=[iindx,i]; iindy=[iindy,j];
else;
iindx=iindx; iindy=iindy;
end;
end;
end;
real_ind=union(iindx,iindy);

%B=G_IB-diag(diag(G_IB)); B=B.*(abs(B)>0.25); i=1; y_d=Y_Data{i};
%B=G_IB-diag(diag(G_IB)); B=B.*(abs(B)>0.15); i=5; y_d=Y_Data{i};

AB2=full(abs(B));  y_d=y_d(:,1)-y_d(:,2);
clf
scatter(diag(AB2(iindx,iindy)), abs(y_d(iindx)-y_d(iindy)))
xlabel('$$\hat{\Omega}_{ij}$$','Interpreter','Latex');
ylabel('$$\Delta_i - \Delta_j$$','Interpreter','Latex');






B=B(real_ind,real_ind);
n=size(B,1);
iindx=[]; iindy=[];
for i=1:n;
for j=1:i;
if abs(B(i,j))>0;
iindx=[iindx,i]; iindy=[iindy,j];
else;
iindx=iindx; iindy=iindy;
end;
end;
end;

clf
G=graph(iindx,iindy);
h = plot(G,'Layout','force');
%labelnode(h, 1:size(real_ind,2),real_ind)
h.NodeLabel=real_ind;
h.MarkerSize=3
set(gca,'fontsize',10)
set(gca,'Visible','off')
%print -depsc Moth_1_IB
print -depsc Moth_5_IB


bins = conncomp(G);
bins2=[real_ind',bins']; 
bin_fin=sortrows(bins2,2); bin_fin=[bin_fin(:,2), bin_fin(:,1)];
bin_num=max(bin_fin(:,1)); bin_siz=[];
for i=1:bin_num
bin_siz=[bin_siz, sum(bin_fin(:,1)==i)];
end
bin_sizz=[(1:bin_num)',bin_siz'];  ord_bin=sortrows(bin_sizz,-2); 

%moth 1
i=1; i=2; i=4
ind=find(bin_fin(:,1)==ord_bin(i,1)); clu_gp=bin_fin(ind,2)
i=3;i=7
i=5; i=9

%moth 5
ind=find(bin_fin(:,1)==ord_bin(i,1)); clu_gp=bin_fin(ind,2); clu_gp'



%moth1
I1=[39,44, 53,107, 112, 117,140, 148,159]; I1=[39 44 112 159]; I1=[107 117 140 148]
I2=[227, 268,283,300,332,389];   I2=[161 227 283]; I2=[268 332 389]
I3=[411, 432, 483, 487,525]; I3=[432 483 487]
X=Total_Data{1};
X1=X(I1,:); X2=X(I2,:); X3=X(I3,:);
YD=Y_Data{1}(:,1)-Y_Data{1}(:,2); [mean(YD(I1)), mean(YD(I2)), mean(YD(I3))], [std(YD(I1)), std(YD(I2)), std(YD(I3))]

%moth 5
I1=[ 12    18    44    69    83   102 ];
I3=[ 208   218   227   237   248   253   263   287];
I2=[104   129   144];
X=Total_Data{5};
X1=X(I1,:); X2=X(I2,:); X3=X(I3,:);
YD=Y_Data{5}(:,1)-Y_Data{5}(:,2); YD=10*YD; [mean(YD(I1)), mean(YD(I2)), mean(YD(I3))], [std(YD(I1)), std(YD(I2)), std(YD(I3))]


clf
h1=plot(X1','b', 'Linewidth', 1); hold on;
h2=plot(X2','k--','Linewidth', 1); hold on
h3=plot(X3','r-.','Linewidth', 1);
xlim([1 500]); ylim([-4 4])
xlabel('Time') %,'FontSize',   23) % x-axis label
ylabel('Torque')
tt=title('Moth L: Wingstrokes');
%print -depsc moth_1_wing
print -depsc moth_5_wing









clf;
figure 
imagesc(G_IB-diag(diag(G_IB))); 
colormap(b2r(-0.1, 0.1)), colorbar;
tt=title('Moth J: $$\hat{\Theta}$$','Interpreter','Latex');
%tt=title(sprintf( ' %s: correlation matrix' ,  title1{i}),'Interpreter','Latex');
set(tt, 'FontSize', 18);
print -depsc Moth_1_IB






