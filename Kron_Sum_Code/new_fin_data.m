
%% this file include moth data anlaysis
load('Total_Dat.mat')    %Total_Datao, Total_Data, Y_Data

title1=cell(1,16);
title1{1}='Moth J'; title1{2}='Moth J';
title1{3}='Moth K'; title1{4}='Moth K';
title1{5}='Moth L'; title1{6}='Moth L';
title1{7}='Moth M'; title1{8}='Moth M';
title1{9}='Moth N'; title1{10}='Moth N';
title1{11}='Moth P'; title1{12}='Moth P';
title1{13}='Moth Q'; title1{14}='Moth Q';
title1{15}='General model'; title1{16}='General model';


%% obtain the temporal inverse covariance matrix estimate \hat{\Theta}

eta=100; R=1;lam=0.01;  

ii=1    %moth J;   ii=5 for moth L
XX=Total_Data{ii}; tauA=1; tauB=1; [n,p]=size(XX);
[G_A, G_IA]=estimator_i2_lassoo(XX,tauA,tauB,eta,lam,R);
%save('Moth_1_graph_A.mat','G_A','G_IA')

%% obtain the spatial inverse covariance matrix estimate \hat{\Omega}

eta=100; R=1;lam=0.1; 
[G_B, G_IB]=estimator_i2_lassoo(XX',tauA,tauB,eta,lam,R);
save('Mooth_1_graph_B.mat','G_B','G_IB')



%% regression
load('Total_Dat.mat')    %Total_Datao, Total_Data, Y_Data

% create a B-spline coefficient matrix
kk=5;
spl_mat=zeros(500,length(-20:3:520)-kk-1);
for x=1:500;
for j=1:(length(-20:3:520)-kk-1)
[vy,vx] = bspline_basis(j,kk,-20:3:520,x); %jntx
spl_mat(x,j)=vy;
end;
end;

spl_mat=spl_mat(:,3:end-2);   %csvwrite('splmat.csv',spl_mat);  
spl_mat=csvread('splmat.csv'); %spl_mat=csvread('splmat.csv');
[U,S,V] = svd(spl_mat,'econ');
CC=V*S*S*V';
C=sqrtm(CC); tauC=trace(C)/size(C,1);
projj=spl_mat*inv(spl_mat'*spl_mat)*spl_mat';


tauA=1;
for ii=1:14;
aaa11=Y_Data{ii}; [n,p]=size(aaa11);
X=Total_Data{ii};  y=aaa11;   y2=y(:,1)-y(:,2); y2=y2-mean(y2); y2=2*sqrt(size(y2,1))*y2/norm(y2);
XX=X*spl_mat; XX=XX-repmat(mean(XX),[n, 1]); X=X-repmat(mean(X),[n, 1]);
tauACB=(norm(XX,'fro')^2)/(size(XX,1)*size(XX,2));
tauA=0.5*tauACB; tauB=(tauACB-tauA)/tauC;
lamr=0.1:3:80; 
parfor jr=1:length(lamr)
b=ridge(y2,X,lamr(jr));   
beta_ridge_d{jr}=b; 
Rs2_ridge(jr)=corr(X*b,y2)^2;
end
bbeta_ridge_d{ii}=beta_ridge_d;
Rsq2_ridge{ii}=Rs2_ridge;

laml=0.01:0.2:4;
parfor jl=1:length(laml)
[B,FitInfo] = lasso(XX,y2,'Lambda', laml(jl)*sqrt(log(size(XX,2))/size(XX,1))); b2=spl_mat*B;
beta_lasso_d{jl}=b2;
Rs2_lasso(jl)=corr(X*b2,y2)^2;
end
bbeta_lasso_d{ii}=beta_lasso_d;
Rsq2_lasso{ii}=Rs2_lasso;
end



load('Moth_1_graph_A.mat');% ii=1
G_A=(G_A+G_A')/2;
test_A=admm(1,G_A, 0.00001, ones(500,500), ones(500,500)); test_A=500*test_A/trace(test_A);

aaa11=Y_Data{ii}; [n,p]=size(aaa11);
X=Total_Data{ii};  y=aaa11;   y2=y(:,1)-y(:,2); y2=y2-mean(y2); y2=2*sqrt(size(y2,1))*y2/norm(y2);
XX=X*spl_mat; XX=XX-repmat(mean(XX),[n, 1]); X=X-repmat(mean(X),[n, 1]);
tauACB=(norm(XX,'fro')^2)/(size(XX,1)*size(XX,2));
tauA=0.5*tauACB; tauB=(tauACB-tauA)/tauC;
[n,p]=size(XX); Gamma=XX'*XX/n-tauB*C;
eta=1*norm(Gamma); R=0.3;  lambda=0.1*sqrt(log(size(XX,2))/size(XX,1)); %ii=1
[est_beta,proc]=composite_regression(tauA,tauB, C, XX,y2, eta,0.1*lambda,R,randn(p,1)); bet=spl_mat*est_beta; 

%%calculate R^2 estimate
denomi=sum((y2-mean(y2)).^2)* n*bet'*test_A*bet
numera=max(y2'*X*bet - norm(y2)*sqrt(n)*norm(bet),0)
EIV_reg_coef_final{ii}=bet;
Rsq2_EIV(ii)=numera^2/denomi








