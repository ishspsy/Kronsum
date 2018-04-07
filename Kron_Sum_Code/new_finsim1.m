%% this file include simulation analysis

%% algorithm requires cvx package in MATLAB
cvx_setup;

%% p=m (in the paper) and d <=p is the sparsity (local maximum sparsity) defined in the paper
p=256; d=2;   

%% generating A=AR(1) matrix using parameter var_AR
    A=zeros(p,p);
    var_AR=0.3;
for i = 1:p ;
    for j = 1:p;
        A(i,j)= var_AR^(abs(i-j));
    end;
end;

%% generating star block model
A=diag(repmat(1,1,p));
AA=repmat(0.25,16,16);   %sub-block matrix
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



%% ss is the set of size of the rescaled sample
ss=[5,9,14,21,30,35,40,45,50,55,60,65,70]; 

%% n_set is the true n dimensions sets
n_set=d^2*log(p)*ss; n_set=round(n_set);

norm22=cell(1,length(n_set));normff=cell(1,length(n_set)); norm22m=zeros(1,length(n_set));normffm=zeros(1,length(n_set));

for ijk=1:length(n_set);
n=n_set(ijk);
mun=zeros(n,1);mup=zeros(p,1);

%% generating B: random model
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

%% trace parameter
tauA=1; tauB=0.3;     % \tau parameter of A and B
A=A*tauA; B=B*tauB;

%% define inverse covariance matrix
IA=inv(A); IB=inv(B);

%% run a simulation

emu=zeros(n,1);
nn=n; pp=p;
norm2=[]; normf=[];

parfor iii=1:100;
Xp=mvnrnd(mup',A,n);
Xn=mvnrnd(mun',B,p);
X=Xp+Xn'; 

% lam1 and R is the regularization parameter defined based on the theory.
% eta is the step-size
R=sqrt(d)*1; eta=1.5*norm(A); 
lam1=0.2*sqrt(log(p)/n)* (sqrt(norm(B))+1); lam2=(sqrt(tauB)+1); lam=lam1*lam2;

%run a nodewise regression% G_A and G_IA estimates A and IA=A^{-1}
[G_A, G_IA]=estimator_i2_lassoo(X,tauA,tauB,eta,lam,R);

%record the performance of the esimate
norm2=[norm2, norm(G_IA-IA)/norm(IA)]; normf=[normf, norm(G_IA-IA,'fro')/norm(IA,'fro')];
end;

norm22{ijk}=norm2; normff{ijk}=normf;
norm22m(ijk)=mean(norm2);
normffm(ijk)=mean(normf) ;
end;

save('fin2_ar03_m5256.mat', 'norm22m', 'normffm' , 'norm22','normff')


