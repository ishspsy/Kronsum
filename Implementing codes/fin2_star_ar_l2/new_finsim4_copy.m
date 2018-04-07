addpath('/home2/sp928/cvx')
addpath('/home2/sp928/MATLAB')
addpath('/home2/sp928/umich_file')

cvx_setup;

%%%%%%%%

%%%%%%512

p=512; d=5;

% AR(1) model
    A=zeros(p,p);
    var_AR=0.3;
for i = 1:p ;
    for j = 1:p;
        A(i,j)= var_AR^(abs(i-j));
    end;
end;


%% star block model
A=diag(repmat(1,1,p));
AA=repmat(0.25,16,16);

for i=2:6;
    AA(1,i)=0.5;
    AA(i,1)=0.5;
end;

for i=0:(round(p/16)-1);
    A((i*16+1):((i+1)*16),(i*16+1):((i+1)*16))=AA;
end;

for i=1:p;
    A(i,i)=1;
end;

ss=[5,9,14,21,30,35,40,45,50,55,60]; n_set=d^2*log(p)*ss; n_set=round(n_set);

norm22=cell(1,length(n_set));normff=cell(1,length(n_set)); norm22m=zeros(1,length(n_set));normffm=zeros(1,length(n_set));

for ijk=1:length(n_set);
n=n_set(ijk);
mun=zeros(n,1);mup=zeros(p,1);

% AR(1) model
    B=zeros(n,n);
    var_AR=0.7;
for i = 1:n ;
    for j = 1:n;
        B(i,j)= var_AR^(abs(i-j));
    end;
end;

ttau_A=1; ttau_B=0.3;  tauA=ttau_A; tauB=ttau_B;
A=A*ttau_A; B=B*ttau_B;

IA=inv(A);
%IB=inv(B);
IAA=(abs(IA)>0.01);
IAA=IAA-diag(diag(IAA));

emu=zeros(n,1);
nn=n; pp=p;
norm2=[]; normf=[];
normB=norm(B);

parfor iii=1:100;

Xp=mvnrnd(mup',A,n);
Xn=mvnrnd(mun',B,p);
X=Xp+Xn'; XX=X;
R=sqrt(d)*1; eta=1.5*norm(A); %lam=0.4*sqrt(log(p)/n);
lam1=0.2*sqrt(log(p)/n)* (sqrt(normB)+1); lam2=(sqrt(tauB)+1); lam=lam1*lam2;
[G_A, G_IA]=estimator_i2_lassoo(XX,tauA,tauB,eta,0.6*lam,R);
norm2=[norm2, norm(G_IA-IA)/norm(IA)]; normf=[normf, norm(G_IA-IA,'fro')/norm(IA,'fro')];
end;
norm22{ijk}=norm2; normff{ijk}=normf;
norm22m(ijk)=mean(norm2);
normffm(ijk)=mean(normf) ;
end;
save('fin2_star_ar07_tauB03_m512.mat', 'norm22m', 'normffm' , 'norm22','normff')

