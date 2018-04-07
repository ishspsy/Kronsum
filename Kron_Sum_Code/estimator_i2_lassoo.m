function[G_A, G_IA]=estimator_i2_lassoo(XX,tauA,tauB,eta,lam,R);

X=XX;
[nn,pp]=size(X);   

%% estimate tauB parameter
tauB=norm(X,'fro')^2/(nn*pp) - tauA;

%% test_A is an initial estimate of A
test_A=XX'*XX/nn-tauB*eye(pp); test_A=pp*tauA*test_A/trace(test_A);
est_A=test_A;

r_est_A=repmat(0,pp,pp); r_est_A1=repmat(0,pp,pp);

for j=1:pp;
y=XX(:,j);
X=XX(:,setdiff(1:pp,j));
[n,p]=size(X);
%ini is the initial vector entering the composite descent algorithm
ini=randn(pp-1,1);
[shat,proc]=composite2_new(tauA,tauB,eye(p),X,y,eta, lam,R,ini);
sibA=ones(pp,1);
%using the relationship of the inverse covariance matrix and the nodewise regression coefficient
sibA(j)=abs(1/(est_A(j,j)-est_A(j,setdiff(1:pp,j))*shat));
sibA(setdiff(1:pp,j))=-abs(1/(est_A(j,j)-est_A(j,setdiff(1:pp,j))*shat))*shat;
r_est_A(:,j)=sibA; 
end;

%% project onto symmetric space by minimizing matrix l1 norm
est_ICorA=proj(r_est_A); 

%% get a estimate of A and A^{-1}
est_CorA=inv(est_ICorA); est_CorA=(est_CorA+est_CorA')/2;
G_A=est_CorA; G_IA=est_ICorA;

