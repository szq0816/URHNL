function [LTrain1,LTrain2,LTrain3,LTrain4] = Noiselabel_function_4dan(LTrain)
% rand('seed',1);
%%
LTrain1 = LTrain;
LTrain2 = LTrain;
LTrain3 = LTrain;
LTrain4 = LTrain;

%%
row = size(LTrain,1);
Rnum4 = round( 0.6 * row);
Rnum3 = round(0.4 * row);
Rnum2 = round(0.2 * row);



%%
Rdata4 = randperm(row,Rnum4); %随机选择ratio*n(n为标签矩阵中1的个数，ratio为选择比率)个数作为缺失标签
Rdata4 = sort(Rdata4);
Rdata3 = Rdata4(1:Rnum3);
Rdata2 = Rdata3(1:Rnum2);




for i = 1:Rnum2
    for j = 1:size(LTrain,2)
        if LTrain(Rdata2(i),j)==1
           LTrain2(Rdata2(i),j) = 0;
           k1 = randperm(size(LTrain,2));
           for k = 1:size(LTrain,2)
               if k1(k)~=j
                   LTrain2(Rdata2(i),k1(k)) = 1;
                   break;
               end 
           end
           break;
        end
    end
end

for i = 1:Rnum3
    for j = 1:size(LTrain,2)
        if LTrain(Rdata3(i),j)==1
           LTrain3(Rdata3(i),j) = 0;
           k1 = randperm(size(LTrain,2));
           for k = 1:size(LTrain,2)
               if k1(k)~=j
                   LTrain3(Rdata3(i),k1(k)) = 1;
                   break;
               end 
           end
           break;
        end
    end
end


for i = 1:Rnum4
    for j = 1:size(LTrain,2)
        if LTrain(Rdata4(i),j)==1
           LTrain4(Rdata4(i),j) = 0;
           k1 = randperm(size(LTrain,2));
           for k = 1:size(LTrain,2)
               if k1(k)~=j
                   LTrain4(Rdata4(i),k1(k)) = 1;
                   break;
               end 
           end
           break;
        end
    end
end


end

