# 数学建模

## 一、评价与决策

### 一、模糊数学模型

评价一个对象优良中差，不能排序

被评判对象的因素论域U，评定等级的论域V，模糊关系矩阵R，评判因素权向量A

A与R合成：B=A*R

模糊算子：

1. M（∧,∨）

```matlab
%函数功能计算两个矩阵的模糊合成
function [R] = synt(A,B)
[m,n] = size(A); [q,p] = size(B);%获得输入的矩阵维度信息
if n ~= q
    disp('第一个矩阵列数和第二个矩阵行数不相同');
else
    R = zeros(m,p); %初始化矩阵
for k = 1:m
    for j = 1:p
        temp = [];
        for i = 1:n
            Min = min(A(k,i),B(i,j));%先取小
            temp = [temp Min];
        end
        R(k,j) = max(temp);%再取大
    end
end
end
end

```

2. M（•,∨）

```matlab
function [R] = synt(A,B)
[m,n] = size(A); [q,p] = size(B);%获得输入的矩阵维度信息
if n ~= q
    disp('第一个矩阵列数和第二个矩阵行数不相同');
else
    R = zeros(m,p); %初始化矩阵
for k = 1:m
    for j = 1:p
        temp = [];
        for i = 1:n
            Mul = A(k,i)*B(i,j);%先相乘
            temp = [temp Mul];
        end
        R(k,j) = max(temp);%再取大
    end
end
end
end
```

3. M（∧,⊕）

```matlab
function [R] = synt(A,B)
[m,n] = size(A); [q,p] = size(B);%获得输入的矩阵维度信息
if n ~= q
    disp('第一个矩阵列数和第二个矩阵行数不相同');
else
    R = zeros(m,p); %初始化矩阵
for k = 1:m
    for j = 1:p
        temp = [];
        for i = 1:n
            Min = min(A(k,i),B(i,j));%先取小
            temp = [temp Min];
        end
        sum = 0;  %再求和
        for i = 1:n
            sum = sum+temp(i);
        end
        R(k,j) = min(1,sum); 
    end
end
end
end

```

4. M（•,⊕）

```matlab
function [R] = synt(A,B)
[m,n] = size(A); [q,p] = size(B);%获得输入的矩阵维度信息
if n ~= q
    disp('第一个矩阵列数和第二个矩阵行数不相同');
else
    R = zeros(m,p); %初始化矩阵
for k = 1:m
    for j = 1:p
        temp = [];
        for i = 1:n
            Mul = A(k,i)*B(i,j); %先相乘
            temp = [temp Mul];
        end
        sum = 0;           %再求和
        for i = 1:n
            sum = sum+temp(i);
        end
        R(k,j) = min(1,sum);
    end
end
end
end

```

## 二、预测

### 一、灰度预测模型——GM(1,1)（尽量不用）

基于原始的数据进行累加计算求得一种规律再进行建模的模型。

优点：将无序的原始序列可以转变为一种有序的生成指数序列。

缺点：在于它只适合于指数增长的预测，较为单一,GM(1,1)为一阶只含一个变量的微分方程模型。 

1. 适用条件

   已知原始样本数量为n ,定义可容覆盖区间Θ

2. GM(1,1)建模公式

   待补充

   [具体参考博客]( https://blog.csdn.net/lvoutongyi/article/details/108037395 )

3. 残差和级别差检验

4. Matlab计算程序

   ```matlab
   function [X0] = HuiDu(X)
   [m,n] = size(X); %获取X矩阵的维度
   X1 = zeros(m,n); %初始化矩阵
   %计算得到X1矩阵
   for j = 1:n
       sum = 0;
       for k = 1:j
           sum = sum + X(k);
       end
       X1(j) = sum;
   end
   
   %计算得到B矩阵
   B = zeros(n-1,2);
   for i = 1:n-1
       B(i,2) = 1;
       B(i,1) = -0.5*(X1(i)+X1(i+1));
   end
   
   %计算得到Y矩阵
   Y = zeros(n-1,1);
   for i = 1:n-1
       Y(i,1) = X(i+1);
   end
   
   %计算得到（a，b）
   A = (inv(B'*B))*B'*Y;
   a = A(1);
   b = A(2);
   
   %计算得到Xc矩阵
   Xc = zeros(1,2*n);
   for i = 0:2*n-1
       Xc(i+1) = (X(1) - b/a) * exp(-a*i) + b/a;
   end
   
   %得到预测矩阵
   X0 = zeros(1,2*n);
   X0(1) = X(1);
   for i = 1:2*n-1
       X0(i+1) = Xc(i+1)-Xc(i); 
   end
   
   end
   
   ```

   

### 二、回归分析模型

### 三、神经网络

### 四、时间序列

1.平稳时间序列模型

[参考博客](https://blog.csdn.net/Will_Zhan/article/details/116425215?ops_request_misc=&request_id=&biz_id=102&utm_term=使用matlab对一组时间序列进行ARMA模型预测&utm_medium=distribute.pc_search_result.none-task-blog-2~all~sobaiduweb~default-3-116425215.first_rank_v2_pc_rank_v29&spm=1018.2226.3001.4187) 

```matlab
%将时间序列平稳化
AZ_diff = zeros(49,5);%必须是列向量
CA_diff = zeros(49,5);
NM_diff = zeros(49,5);
TX_diff = zeros(49,5);
%得到了差分之后的向量
for i = 1:5
    AZ_diff(:,i) = (diff(AZ(i,:)))';%转置为列向量
    CA_diff(:,i) = (diff(CA(i,:)))';%转置为列向量
    NM_diff(:,i) = (diff(NM(i,:)))';%转置为列向量
    TX_diff(:,i) = (diff(TX(i,:)))';%转置为列向量
end

%时间序列的平稳化与ARMA模型估计
% ARMA_model
%估计阶数建立模型
AR_Order = 2;
MA_Order = 2;
MA1 = arima(AR_Order, 0, MA_Order);
%获取参数
step = 41;%需要预测的个数。
%AZ
auto_fore_AZ = zeros(step,5);%初始化
for i = 1:5
    EstMdl = estimate(MA1,AZ_diff(:,i));
    auto_fore_AZ(:,i) = forecast(EstMdl,step,'Y0',AZ_diff(:,i));%自动估计今后step年的数据
end
%CA
auto_fore_CA = zeros(step,5);%初始化
for i = 1:5
    EstMdl = estimate(MA1,CA_diff(:,i));
    auto_fore_CA(:,i) = forecast(EstMdl,step,'Y0',CA_diff(:,i));%自动估计今后step年的数据
end
%NM
auto_fore_NM = zeros(step,5);%初始化
for i = 1:5
    EstMdl = estimate(MA1,NM_diff(:,i));
    auto_fore_NM(:,i) = forecast(EstMdl,step,'Y0',NM_diff(:,i));%自动估计今后step年的数据
end
%TX
auto_fore_TX = zeros(step,5);%初始化
for i = 1:5
    EstMdl = estimate(MA1,TX_diff(:,i));
    auto_fore_TX(:,i) = forecast(EstMdl,step,'Y0',TX_diff(:,i));%自动估计今后step年的数据
end

%得到预测的实际能量值
AZ_fore = zeros(step,5);
CA_fore = zeros(step,5);
NM_fore = zeros(step,5);
TX_fore = zeros(step,5);
%AZ
for i = 1:5
    AZ_fore(1,i) = AZ(i,50) + auto_fore_AZ(1,i); %2010年的数据
    for j = 2:step
        AZ_fore(j,i) = AZ_fore(j-1,i) + auto_fore_AZ(j,i);
    end
end
%CA
for i = 1:5
    CA_fore(1,i) = CA(i,50) + auto_fore_CA(1,i); %2010年的数据
    for j = 2:step
        CA_fore(j,i) = CA_fore(j-1,i) + auto_fore_CA(j,i);
    end
end
%NM
for i = 1:5
    NM_fore(1,i) = NM(i,50) + auto_fore_NM(1,i); %2010年的数据
    for j = 2:step
        NM_fore(j,i) = NM_fore(j-1,i) + auto_fore_NM(j,i);
    end
end
%TX
for i = 1:5
    TX_fore(1,i) = TX(i,50) + auto_fore_TX(1,i); %2010年的数据
    for j = 2:step
        TX_fore(j,i) = TX_fore(j-1,i) + auto_fore_TX(j,i);
    end
end
```







## 参考文章

[如何根据实际问题选择一个合适的数学模型](https://blog.csdn.net/weixin_43014927/article/details/99701935) 















