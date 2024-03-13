function [map1, map2] = URHNL(X1, X2, Y1, Y2, param, A, X1_UR, X1_UQ, X2_UR, X2_UQ, L_UR, L_UQ, LX1, LY1 )

alphe1 = param.alphe1; alphe2 = param.alphe2;
beta1 = param.beta1; beta2 = param.beta2;
lambda = param.lambda; thea = param.thea;
gamma = param.gamma; mu = param.mu;
lambda1 = param.lambda1; alphe3 = param.alphe3;
mu1 = param.mu1; lambda_c = param.lambda_c;


% tic;

X1 = X1'; X2=X2'; Y1 = Y1'; Y2 = Y2';
[dx1, n ] = size(X1); [dx2, ~ ] = size(X2); [c, ~] = size(Y1);


rho = 1;
etaX = 0.5; etaY = 0.5; 

a = (c*(c+2)+c*sqrt(c*(c+2)))/4+eps; %More precise value can be experimentally obtained.
r = param.nbits;



%% 初始化X1的噪声标签相关参数
D1 = zeros(size(Y1));
Yn1 = Y1;
F11 = zeros(size(Y1));
K1 = zeros(size(Y1));
Y11 = zeros(size(Y1));

Y01 = Y1;
ind11 = Y01>0;
ind12 = Y01 == 0;
k = inf;


%% 初始化X1的噪声标签相关参数
D2 = zeros(size(Y2));
Yn2 = Y2;
F12 = zeros(size(Y2));
K2 = zeros(size(Y2));
Y12 = zeros(size(Y2));

Y02 = Y2;
ind21 = Y02>0;
ind22 = Y02 == 0;

%%
sel_sample = X1(:,randsample(n, 500),:);
[pcaW, ~] = eigs(cov(sel_sample'), r);
V1_c = pcaW'*X1;
sel_sample = X2(:,randsample(n, 500),:);
 [pcaW, ~] = eigs(cov(sel_sample'), r);
V2_c = pcaW'*X2;
V1 =  (V1_c+V2_c)/2 ;
V2 =  (V1_c+V2_c)/2 ;




t1 = ones(r,1);
t2 = ones(r,1);
en = ones(1,n);
P1 = rand(r,c);
P2 = rand(r,c);


%%
AA = A'*A;
YY1 = Yn1*Yn1';
YY2 = Yn2*Yn2';
X1X1 = X1*X1';
X2X2 = X2*X2'; 



 for iter = 1:param.iter 
     
      %% B1 and B2
      B1 = sign(lambda*V1 + (r*thea)*(V1*Yn1'*Yn1));
      B1(B1==0)=-1;
      
      B2 = sign(lambda*V2 +(r*thea)*(V2*Yn2'*Yn2));
      B2(B2==0)=-1;
     
     %% Z1 and Z2
      Z1 = pinv(alphe1*AA)*(alphe1*A'*X1*Yn1')*pinv(YY1);
      Z2 = pinv(alphe2*AA)*(alphe1*A'*X2*Yn2')*pinv(YY2);
      
      %% H1 and H2
      % 正交约束
      [U11 ,~ ,V111] = svd(alphe3*V1*Yn1',0);
       H1= U11  * V111';

       [U22 ,~ ,V222] = svd(alphe3*V2*Yn2',0);
       H2= U22  * V222';
      
      
     %% Yn1 and Yn2
      Yn1 = pinv(alphe1*Z1'*A'*A*Z1+(alphe3+mu1)*eye(c))*(rho/2*(Y1-D1+1/rho*F11)+alphe1*Z1'*A'*X1+alphe3*H1'*V1+(mu1*K1-Y11));
      Yn1(ind11) = e_dragging(Yn1(ind11),1,k);
      Yn1(ind12) = e_dragging(Yn1(ind12),-k,0);
      
      Yn2 = pinv(alphe2*Z2'*A'*A*Z2+(alphe3+mu1)*eye(c))*(rho/2*(Y2-D2+1/rho*F12)+alphe2*Z2'*A'*X2+alphe3*H2'*V2+(mu1*K2-Y12));
      Yn2(ind21) = e_dragging(Yn2(ind21),1,k);
      Yn2(ind22) = e_dragging(Yn2(ind22),-k,0);
      
      %% K1 and K2
      [U11_1,S11,V11_1] = svd(Yn1 + Y11./mu1,'econ');
      a11 = diag(S11)-lambda_c/mu1;
      a11(a11<0)=0; 
      T1 = diag(a11);
      K1 = U11_1*T1*V11_1'; 
      
      [U12,S12,V11_2] = svd(Yn2 + Y12./mu1,'econ');
      a12 = diag(S12)-lambda_c/mu1;
      a12(a12<0)=0; 
      T2 = diag(a12);
      K2 = U12*T2*V11_2'; 
     
      %% D1 and D2
      Etp1 = Y1 - Yn1 + 1/rho*F11;
      D1 = sign(Etp1).*max(abs(Etp1)- mu/rho,0); 
      
      Etp2 = Y2 - Yn2 + 1/rho*F12;
      D2 = sign(Etp2).*max(abs(Etp2)- mu/rho,0); 
      
      %% F11 and F12
      F11 = F11 + rho*(Y1 - Yn1 -D1);
      
      F12 = F12 + rho*(Y2 - Yn2 -D2);
      rho = min(1e4, 1.3*rho);
     
      %% P1 and P2
      P1 = (beta1*V1-beta1*t1*en)*X1'*pinv(beta1*X1X1+gamma*eye(dx1));
      P2 = (beta2*V2-beta2*t2*en)*X2'*pinv(beta2*X2X2+gamma*eye(dx2));
      
      %% t1 and t2
      t1 = (1/n)*(V1-P1 *X1)*en';
      t2 = (1/n)*(V2-P2 *X2)*en';
      
      
      %% V
      V_1 = beta1*(P1*X1+t1*en)+(r*thea)*(B1*Yn1'*Yn1)+ (r*lambda1)/(a+etaX)*(a*V2*Yn1'*Yn2+etaX*V2*X1'*X2)+ lambda*B1 + alphe3*H1*Yn1;
      V1 = solve_V(V_1,n,r);
       
      V_2 = beta2*(P2*X2+t2*en)+(r*thea)*(B2*Yn2'*Yn2)+(r*lambda1)/(a+etaX)*(a*V1*Yn2'*Yn1+etaX*V1*X2'*X1)+lambda*B2 + alphe3*H2*Yn2;
      V2 = solve_V(V_2,n,r);
      
      
      %% Y11 and Y12
      mu1 = 1.01*mu1;
      Y11 = Y11 + mu1*(Yn1-K1);
      
      Y12 = Y12 + mu1*(Yn2-K2);
      

%% Compute mAP
    rBX = [sign(P1 * X1_UR' - t1 * ones(1,size(X1_UR,1)))';B1'];
    qBX = sign(P1 * X1_UQ' - t1 * ones(1,size(X1_UQ,1)))';
    rBY = [sign(P2 * X2_UR' - t2 * ones(1,size(X2_UR,1)))';B2'];
    qBY = sign(P2 * X2_UQ' - t2 * ones(1,size(X2_UQ,1)))';

    rBX = (rBX > 0);
    qBX = (qBX > 0);
    rBY = (rBY > 0);
    qBY = (qBY > 0);
    
    fprintf('\nText-to-Image Result:\n');
    mapTI = map_rank1(qBY, rBX,[L_UR; LX1], L_UQ);
    map1(iter) = mapTI(end);

    
    fprintf('Image-to-Text Result:\n');
    mapIT = map_rank1(qBX, rBY,[L_UR; LY1], L_UQ);
    map2(iter) =  mapIT(end);
      
 end
 end

