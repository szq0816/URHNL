clear;clc

addpath data;
addpath utils;

 %% 数据处理
 
  %% load data
    load 'Labelme_label.mat'
    load 'Labelme_XY.mat'
    load 'Labelme_label_noise_X.mat'
    load 'Labelme_label_noise_Y.mat'
    load 'Labelme_attributes.mat'

    LX1 = LX_SR;
    LY1 = LY_SR;
    S = labelme_attributes;
    
    
run = 10;
for i = 1 : run

%% RZSDH 不配对参数
param.alphe1 = 1e1; param.alphe2 = 1e1;
param.beta1  = 1e4; param.beta2  = 1e4;
param.thea  = 1e-3; param.lambda = 1e-3;
param.mu = 1e-3; param.gamma = 1e2;

 
param.lambda_c = 1e-1; param.mu1 = 1e-4;
param.lambda1 = 1e-3; param.alphe3 = 1e-4;


param.iter = 10;
nbitset  = 32;

eva_info = cell(1,length(nbitset));
for bit = 1:length(nbitset) 
    %% 无噪声
    param.nbits = nbitset(bit);
    [map11, map22] = URHNL(X, Y, LTrain11, LTrain21, param, S', X1_UR, X1_UQ, X2_UR, X2_UQ, L_UR, L_UQ, LX1, LY1 );

    map1(i,1) = max(map11);
    map1(i,2) = max(map22);
   
   

 %% 0.2的噪声
    param.nbits = nbitset(bit);
    [map11, map22] = URHNL(X, Y, LTrain12, LTrain22, param, S', X1_UR, X1_UQ, X2_UR, X2_UQ, L_UR, L_UQ, LX1, LY1 );

    map2(i,1) = max(map11);
    map2(i,2) = max(map22);


 %% 0.4的噪声
    param.nbits = nbitset(bit);
    [map11, map22] = URHNL(X, Y, LTrain13, LTrain23, param, S', X1_UR, X1_UQ, X2_UR, X2_UQ, L_UR, L_UQ, LX1, LY1 );

    map3(i,1) = max(map11);
    map3(i,2) = max(map22);
    


 %% 0.6的噪声
    param.nbits = nbitset(bit);
    [map11, map22] = URHNL(X, Y, LTrain14, LTrain24, param, S', X1_UR, X1_UQ, X2_UR, X2_UQ, L_UR, L_UQ, LX1, LY1 );


    map4(i,1) = max(map11);
    map4(i,2) = max(map22);

end
end
map111 = mean(map1)
map222 = mean(map2)
map333 = mean(map3)
map444 = mean(map4)

