addpath('/home2/sp928/cvx')
addpath('/home2/sp928/MATLAB')
addpath('/home2/sp928/umich_file')

cvx_setup;

p=128; d=2;

% AR(1) model
    A=zeros(p,p);
    var_AR=0.3;
for i = 1:p ;
    for j = 1:p;
        A(i,j)= var_AR^(abs(i-j));
    end;
end;

ss=[5,9,14,21,30,35,40,45,50,55,60,65,70]; n_set=d^2*log(p)*ss; n_set=round(n_set);
norm22=cell(1,length(n_set));normff=cell(1,length(n_set)); norm22m=zeros(1,length(n_set));normffm=zeros(1,length(n_set));

for ijk=1:length(n_set);
n=n_set(ijk);
mun=zeros(n,1);mup=zeros(p,1);

%B: random model
set= [];
for k=0:(n-2);
    set=[set,(k*n+(k+2)):((k+1)*n)];
end;
chose=datasample(set,round(1*n),'Replace',false);
B=diag(repmat(0.25,1,n));
for k=1:round(n*1);
    rand=unifrnd(0.1,0.3);
    ind1=floor((chose(k)-1)/n);
    ind2=chose(k)-ind1*n;
    ind1=ind1+1;
    B(ind2,ind1)=B(ind2,ind1)-rand;
    B(ind1,ind2)=B(ind1,ind2)-rand;
    B(ind1,ind1)=B(ind1,ind1)+rand;
  B(ind2,ind2)=B(ind2,ind2)+rand;
end;

B=inv(B);
B=B*n/trace(B);

ttau_A=1; ttau_B=0.5;  tauA=ttau_A; tauB=ttau_B;
A=A*ttau_A; B=B*ttau_B;

IA=inv(A);
IB=inv(B);
IAA=(abs(IA)>0.01);
IAA=IAA-diag(diag(IAA));

emu=zeros(n,1);
nn=n; pp=p;
norm2=[]; normf=[];

nrep=4;
parfor iii=1:100;

Xp=mvnrnd(mup',A,n); XXN=cell(1,nrep);
for df=1:nrep;
XXN{df}=Xp+mvnrnd(mun',B,p)';
end;
XXX=0;
for df=1:nrep;
XXX=XXX+XXN{df}/nrep;
end
XX=XXX; X=XX;
R=sqrt(d)*1; eta=1.5*norm(A); %lam=0.4*sqrt(log(p)/n);
lam1=0.2*sqrt(log(p)/n)* (sqrt(norm(B/nrep))+1); lam2=(sqrt(tauB/nrep)+1); lam=lam1*lam2;
[G_A, G_IA]=estimator_i2_lassoo(XX,tauA,tauB/nrep,eta,0.6*lam,R);
norm2=[norm2, norm(G_IA-IA)/norm(IA)]; normf=[normf, norm(G_IA-IA,'fro')/norm(IA,'fro')];
end;
norm22{ijk}=norm2; normff{ijk}=normf;
norm22m(ijk)=mean(norm2);
normffm(ijk)=mean(normf) ;
end;

save('fin2_ar03_tauB05_N4_m128.mat', 'norm22m', 'normffm' , 'norm22','normff')

%save('fin_ar03_N2_m128.mat', 'norm22m', 'normffm' , 'norm22','normff')


%%%%%%%%

p=256; d=2;

% AR(1) model
    A=zeros(p,p);
    var_AR=0.3;
for i = 1:p ;
    for j = 1:p;
        A(i,j)= var_AR^(abs(i-j));
    end;
end;

ss=[5,9,14,21,30,35,40,45,50,55,60,65,70]; n_set=d^2*log(p)*ss; n_set=round(n_set);
norm22=cell(1,length(n_set));normff=cell(1,length(n_set)); norm22m=zeros(1,length(n_set));normffm=zeros(1,length(n_set));

for ijk=1:length(n_set);
n=n_set(ijk);
mun=zeros(n,1);mup=zeros(p,1);

%B: random model
set= [];
for k=0:(n-2);
    set=[set,(k*n+(k+2)):((k+1)*n)];
end;
chose=datasample(set,round(1*n),'Replace',false);
B=diag(repmat(0.25,1,n));
for k=1:round(n*1);
    rand=unifrnd(0.1,0.3);
    ind1=floor((chose(k)-1)/n);
    ind2=chose(k)-ind1*n;
    ind1=ind1+1;
    B(ind2,ind1)=B(ind2,ind1)-rand;
    B(ind1,ind2)=B(ind1,ind2)-rand;
    B(ind1,ind1)=B(ind1,ind1)+rand;
  B(ind2,ind2)=B(ind2,ind2)+rand;
end;

B=inv(B);

B=B*n/trace(B);

ttau_A=1; ttau_B=0.5;  tauA=ttau_A; tauB=ttau_B;
A=A*ttau_A; B=B*ttau_B;

IA=inv(A);
IB=inv(B);
IAA=(abs(IA)>0.01);
IAA=IAA-diag(diag(IAA));

emu=zeros(n,1);
nn=n; pp=p;
norm2=[]; normf=[];

nrep=4;
parfor iii=1:100;

Xp=mvnrnd(mup',A,n); XXN=cell(1,nrep);
for df=1:nrep;
XXN{df}=Xp+mvnrnd(mun',B,p)';
end;
XXX=0;
for df=1:nrep;
XXX=XXX+XXN{df}/nrep;
end
XX=XXX;	X=XX;
R=sqrt(d)*1; eta=1.5*norm(A); %lam=0.4*sqrt(log(p)/n);
lam1=0.2*sqrt(log(p)/n)* (sqrt(norm(B/nrep))+1); lam2=(sqrt(tauB/nrep)+1); lam=lam1*lam2;
[G_A, G_IA]=estimator_i2_lassoo(XX,tauA,tauB/nrep,eta,0.6*lam,R);
norm2=[norm2, norm(G_IA-IA)/norm(IA)]; normf=[normf, norm(G_IA-IA,'fro')/norm(IA,'fro')];
end;
norm22{ijk}=norm2; normff{ijk}=normf;
norm22m(ijk)=mean(norm2);
normffm(ijk)=mean(normf) ;
end;

save('fin2_ar03_tauB05_N4_m256.mat', 'norm22m', 'normffm' , 'norm22','normff')

%save('fin_ar03_N2_m256.mat', 'norm22m', 'normffm' , 'norm22','normff')


%%%%%%512

p=512; d=2;

% AR(1) model
    A=zeros(p,p);
    var_AR=0.3;
for i = 1:p ;
    for j = 1:p;
        A(i,j)= var_AR^(abs(i-j));
    end;
end;

ss=[5,9,14,21,30,35,40,45,50,55,60,65,70]; n_set=d^2*log(p)*ss; n_set=round(n_set);
norm22=cell(1,length(n_set));normff=cell(1,length(n_set)); norm22m=zeros(1,length(n_set));normffm=zeros(1,length(n_set));

for ijk=1:length(n_set);
n=n_set(ijk);
mun=zeros(n,1);mup=zeros(p,1);

%B: random model
set= [];
for k=0:(n-2);
    set=[set,(k*n+(k+2)):((k+1)*n)];
end;
chose=datasample(set,round(1*n),'Replace',false);
B=diag(repmat(0.25,1,n));
for k=1:round(n*1);
    rand=unifrnd(0.1,0.3);
    ind1=floor((chose(k)-1)/n);
    ind2=chose(k)-ind1*n;
    ind1=ind1+1;
    B(ind2,ind1)=B(ind2,ind1)-rand;
    B(ind1,ind2)=B(ind1,ind2)-rand;
    B(ind1,ind1)=B(ind1,ind1)+rand;
  B(ind2,ind2)=B(ind2,ind2)+rand;
end;

B=inv(B);

B=B*n/trace(B);

ttau_A=1; ttau_B=0.5;  tauA=ttau_A; tauB=ttau_B;
A=A*ttau_A; B=B*ttau_B;

IA=inv(A);
IB=inv(B);
IAA=(abs(IA)>0.01);
IAA=IAA-diag(diag(IAA));

emu=zeros(n,1);
nn=n; pp=p;
norm2=[]; normf=[];

nrep=4;
parfor iii=1:100;

Xp=mvnrnd(mup',A,n); XXN=cell(1,nrep);
for df=1:nrep;
XXN{df}=Xp+mvnrnd(mun',B,p)';
end;
XXX=0;
for df=1:nrep;
XXX=XXX+XXN{df}/nrep;
end
XX=XXX; X=XX;
R=sqrt(d)*1; eta=1.5*norm(A); %lam=0.4*sqrt(log(p)/n);
lam1=0.2*sqrt(log(p)/n)* (sqrt(norm(B/nrep))+1); lam2=(sqrt(tauB/nrep)+1); lam=lam1*lam2;
[G_A, G_IA]=estimator_i2_lassoo(XX,tauA,tauB/nrep,eta,0.6*lam,R);
norm2=[norm2, norm(G_IA-IA)/norm(IA)]; normf=[normf, norm(G_IA-IA,'fro')/norm(IA,'fro')];
end;
norm22{ijk}=norm2; normff{ijk}=normf;
norm22m(ijk)=mean(norm2);
normffm(ijk)=mean(normf) ;
end;

save('fin2_ar03_tauB05_N4_m512.mat', 'norm22m', 'normffm' , 'norm22','normff')


%save('fin_ar03_N2_m512.mat', 'norm22m', 'normffm' , 'norm22','normff')

