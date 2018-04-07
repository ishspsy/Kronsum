%% gradient descent algorithm
function[est_beta,proc]=composite2_new(tauA,tauB, C, X,y,eta,lambda,R,b);

%% We consider EIV regression
% y=X0 beta + epsilon, where X=X0+W, where vec(X0)~N(0, A \otimes I) and vec(W)~N(C \otimes B) by solving 
%%% min_{beta} 1/2 beta' Gamma beta - <r,beta> + lambda |beta|_1
%%% s.t |beta|_1 <= R

%% input 
% X: design matrix
% y: respose vector
% eta: step-size
% lambda: penalty parameter
% R: penalty parameter imposed on the constraint
% b: initial input in the algorithm

[n,p]=size(X);  %in theory, we use notations f=n, m=p

%hat_B is \hat{{tr}}(B)/n
hat_B=tauB;

%Gamma
Gamma=X'*X/n -C*hat_B;

%r
r=X'*y/n;

%initial setting
bb=ones(p,1);

proc=cell(1,3000);
iter=0;

%% repeat until converge
while((norm(bb-b,1)>0.001)+(iter<3000)>1.5);
bb=b;
iter=iter+1;
v=bb-(Gamma*bb-r)/eta;
v=times(v-lambda,(v>lambda))+times(v+lambda,(v<-lambda)); %soft-thresholding
b=projl1(v,R);
proc{iter}=b;
end;

est_beta=b;
