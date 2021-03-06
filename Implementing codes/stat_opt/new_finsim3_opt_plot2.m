addpath('/home2/sp928/cvx')
addpath('/home2/sp928/MATLAB')
addpath('/home2/sp928/umich_file')

cvx_setup;

% STAR: p=512; n=256;  d=5;

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

Xp=mvnrnd(mup',A,n);
Xn=mvnrnd(mun',B,p);
X=Xp+Xn'; XX=X;
R=sqrt(d)*1; eta=1.5*norm(A); %lam=0.4*sqrt(log(p)/n);
lam1=0.2*sqrt(log(p)/n)* (sqrt(normB)+1); lam2=(sqrt(tauB)+1); lam=lam1*lam2;
[nn, pp]=size(XX);

X=XX;

tauB=norm(X,'fro')^2/(nn*pp) - tauA;
test_A=XX'*XX/nn-tauB*eye(pp); test_A=pp*tauA*test_A/trace(test_A);
est_A=test_A;

%ran_set=sort(randperm(pp,20));

ran_set1=1:16:pp;    %herb
  %ran_set2=sort(randperm(setdiff(pp,ran_set1),20));  %leaf
ran_set2=[2,3,4,5,6];ran_set2=[ran_set2,ran_set2+16,ran_set2+32,ran_set2+48];            %leaf


ran_set=ran_set1;  %%herb
ran_set=ran_set2;  %%leaf  

clear var sstat1; clear var sstat2; clear var oopt1; clear var oopt2;

for jk=1:10 
iind=0;

for j=ran_set
iind=iind+1;
set_di=setdiff(1:pp,j);
theta_j=inv(A(set_di,set_di))*A(set_di,j);

y=XX(:,j);
X=XX(:,setdiff(1:pp,j));
[n,p]=size(X);

ini=randn(pp-1,1);
[shat,proc]=composite2_inf(tauA,tauB,eye(p),X,y,eta, lam,R,ini);

sstat=[]; oopt=[];
for ii=1:100
sstat=[sstat,norm(proc{ii}-theta_j)];
oopt=[oopt,norm(proc{ii}-shat)];
end;
sstat=log(sstat); oopt=log(oopt);

sstat2{iind}=sstat; oopt2{iind}=oopt;
end

sstat1{jk}=sstat2; oopt1{jk}=oopt2;
end


sstat3=cell(1,length(ran_set)); oopt3=cell(1,length(ran_set));
for ii=1:length(ran_set)

st3=0; op3=0; 
for jk=1:10
st3=st3+sstat1{jk}{ii}/10;
op3=op3+oopt1{jk}{ii}/10;
end
sstat3{ii}=st3;
oopt3{ii}=op3;
end

% sstat3_herb=sstat3; oopt3_herb=oopt3;
% sstat3_leaf=sstat3; oopt3_leaf=oopt3;


save('opt_stat_STAR.mat', 'sstat3_herb', 'oopt3_herb','sstat3_leaf', 'oopt3_leaf')
