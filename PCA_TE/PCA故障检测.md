（1）     计算T2统计量控制限和SPE统计量控制限；

​	![1572071054425](C:\Users\C&Y0903\AppData\Roaming\Typora\typora-user-images\1572071054425.png)



（2 ）计算测试样本的T2统计量和SPE统计量

![1572071114734](C:\Users\C&Y0903\AppData\Roaming\Typora\typora-user-images\1572071114734.png)



```matlab
data=load('TE_data.mat');
%数据读取
data = struct2cell(data);
testdata =  data(1:22);
train = cell2mat(data(23));
train = train';
train_mean = mean(train);  %按列 Xtrain 平均值                           
train_std = std(train);    %求标准差                       
[train_row,train_col] = size(train); %求 train 行、列数                                                                
train=(train-repmat(train_mean,train_row,1))./repmat(train_std,train_row,1); 
%归一化
%求协方差矩阵 
sigmatrain = cov(train);
%对协方差矩阵进行特征分解，lamda 为特征值构成的对角阵，T的列为单位特征向量，且与 lamda 中的特征值一一对应：
[T,lamda] = eig(sigmatrain);
disp('特征根（由小到大）'); 
disp(lamda); 
disp('特征向量：'); 
disp(T);
%取对角元素(结果为一列向量)，即 lamda 值，并上下反转使其从大到小排列，主元个数初值为 1，若累计贡献率小于 90
%则增加主元个数 
D = flipud(diag(lamda));                             
num_pc = 1;                                          
while sum(D(1:num_pc))/sum(D) < 0.9    
num_pc = num_pc +1; 
end   
%取与 lamda 相对应的特征向量 
P = T(:,train_col-num_pc+1:train_col);
%每一列代表一个特征向量
%求置信度为 99%、95%时的 T2 统计控制限                        
T2UCL1=num_pc*(train_row-1)*(train_row+1)*finv(0.99,num_pc,train_row - num_pc)/(train_row*(train_row - num_pc)); 
T2UCL2=num_pc*(train_row-1)*(train_row+1)*finv(0.95,num_pc,train_row - num_pc)/(train_row*(train_row - num_pc)); 
%开始计算SPE统计量
for i = 1:3 
    theta(i) = sum((D(num_pc+1:train_col)).^i); 
end 
h0 = 1 - 2*theta(1)*theta(3)/(3*theta(2)^2); 
ca = norminv(0.99,0,1); 
SPE = theta(1)*(h0*ca*sqrt(2*theta(2))/theta(1) + 1 + theta(2)*h0*(h0 - 1)/theta(1)^2)^(1/h0); 
%计算出了SPE统计量的界限
for k = 1:22
    %22组测试数据
    test = cell2mat(testdata(k));
    %开始在线检测
    n = size(test,1); 
    test=(test-repmat(train_mean,n,1))./repmat(train_std,n,1); 
    %测试样本归一化
    [r,y] = size(P*P'); 
    I = eye(r,y); %单位矩阵
    T2_test = zeros(n,1); 
    SPE_test = zeros(n,1); 
    for i = 1:n 
        T2_test(i)=test(i,:)*P*inv(lamda(52-num_pc+1:52,52-num_pc+1:52))*P'*test(i,:)';                                           
        SPE_test(i) = test(i,:)*(I - P*P')*test(i,:)';                                                                                    
    end 
    %绘图 
    figure (k);
    subplot(2,1,1); 
    plot(1:n,T2_test,'k');                                     
    title('主元分析统计量变化图T2'); 
    xlabel('采样数'); 
    ylabel('T^2'); 
    hold on;     
    line([0,n],[T2UCL1,T2UCL1],'LineStyle','--','Color','r');%画出标志线 
    line([0,n],[T2UCL2,T2UCL2],'LineStyle','--','Color','g'); 
    subplot(2,1,2); 
    plot(1:n,SPE_test,'k'); 
    title('主元分析统计量变化图SPE')
    xlabel('采样数'); 
    ylabel('SPE'); 
    hold on;  
    line([0,n],[SPE,SPE],'LineStyle','--','Color','r'); 
end


```

仿真结果图：

![1572070843447](C:\Users\C&Y0903\AppData\Roaming\Typora\typora-user-images\1572070843447.png)

这是测试te_00的仿真结果。

上述代码运行会输出所有数据的仿真结果，（这里粘贴了te_00为例子）

