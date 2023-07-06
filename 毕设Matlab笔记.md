---
title: 毕设Matlab笔记
tags: [Matlab]
categories: Matlab
index_img: /img/default.png
banner_img: /img/default.png
---

毕设Matlab笔记

<!-- more -->

# Matlab 笔记

## 一、PA建模

本例使用NXP Airfast LDMOS Doherty PA，工作频率3.6-3.8 GHz，增益29 dB。

## 1、无记忆非线性模型阅读

导入原始测量的功放输入输出数据（复数形式）

对输入输出数据取绝对值，然后转换为dBm

对数据进行统计histcounts，直方图

求输入数据最大值-20

找出满足条件的边缘数据的位置

建立新的表格存储满足条件的数据值

分别求输入功率，输出功率，相移的平均值

```matlab
sampleRate = 860160000;

%数据转换
paInputdBm  = mag2db(abs(paInput)) + 30 - 20;
paOutputdBm  = mag2db(abs(paOutput)) + 30 - 20;

%准备建模数据
[N,edges,idx] = histcounts(paInputdBm, 'BinWidth', 0.5);  %hiscounts函数，直方图bin计数
minInPowerdBm = max(paInputdBm) - 20;  %输入功率值不能小于最大输入功率20db
minIdx = find(edges < minInPowerdBm, 1, 'last');
tableLen = length(edges)-minIdx-1;
inOutTable = zeros(tableLen,2);
for p = minIdx+1:length(edges)-1
	inOutTable(p-minIdx,1) = mean(paInputdBm(idx == p));   % Average input power for current bin
	inOutTable(p-minIdx,2) = mean(paOutputdBm(idx == p));  % Average output power for current bin
	inOutTable(p-minIdx,3) = mean(angle(paOutput(idx == p)./paInput(idx == p))); % Average phase shift for current                                                                                       bin
end

%建立模型
pa = comm.MemorylessNonlinearity('Method','Lookup table','Table',inOutTable,'ReferenceImpedance',100);

%计算模型输出
paOutputFitMemless = pa(paInput);
%计算误差
err = abs(paOutput - paOutputFitMemless)./abs(paOutput);
rmsErrorMemless = rms(err)*100;
disp(['Percent RMS error in time domain is ' num2str(rmsErrorMemless) '%']);
%绘制图形观察
helperPACharPlotTime(paOutput, paOutputFitMemless, sampleRate);
helperPACharPlotGain(paInput, paOutput, paOutputFitMemless);

```



## 2、记忆多项式模型阅读



程序阅读

```matlab
modType = 'memPoly'; %记忆多项式模型
memLen = 5;     % M = 5
degLen = 5;     % K = 5
numDataPts = length(paInput);
halfDataPts = round(numDataPts/2);
fitCoefMatMem = helperPACharMemPolyModel('coefficientFinder',paInput(1:halfDataPts),paOutput(1:halfDataPts),memLen,degLen,modType);
disp(abs(fitCoefMatMem));

rmsErrorTimeMem = helperPACharMemPolyModel('errorMeasure',paInput, paOutput, fitCoefMatMem, modType);
disp(['Percent RMS error in time domain is ' num2str(rmsErrorTimeMem) '%']);

paOutputFitMem = helperPACharMemPolyModel('signalGenerator', paInput, fitCoefMatMem, modType);

helperPACharPlotTime(paOutput, paOutputFitMem, sampleRate);
helperPACharPlotGain(paInput, paOutput, paOutputFitMem);



```

我的程序

```matlab
%x = (0 : 3);
x = paInput;
x = x(:);
xLength = length(x);

%y = (2 : 5);
y = paOutput;
y = y(:);
yLength = length(y);

M = 3;
K = 5;

%构造(x(n) - m)矩阵，m从0到M-1, 矩阵大小为xLength*M
xm = x;
for m = 0 : M - 1
    xm(M:xLength , m+1) = xm(M-m:xLength-m , 1) ;
end
xm_fix = xm; %保存(x(n) - m)矩阵

%构造(x(n) - m) * |x(n) - m|^k 矩阵，矩阵大小为(xLength-M+1)*（M*K）
xmAbs = abs(xm);%保存|x - m|矩阵
for k = 1 : K-1
    mid = (xmAbs.^ k) .* xm_fix;
    xm = cat(2, xm, mid);
end
xmk = xm(M:xLength,:);

%计算系数矩阵，矩阵大小为（M*K） * 1
coef = xmk \ y(M:xLength);

%带入求解yout矩阵，矩阵大小为xLength*1
yout = xm * coef;

% helperPACharPlotTime(y, yout, 860160000);
% helperPACharPlotGain(x, y, yout);


%通过归一化均方误差衡量功放的建模精度
NMSE = 10 * log10( sum ( (abs(y - yout)).^2 ) / sum ( (abs(y).^2) ) ) ;



%AM/AM图绘制
%单位转换
paInputPowerdBm = mag2db(abs(x)) + 30 - 20;
paOutputPowerdBm = mag2db(abs(y)) + 30 - 20;
paOutputPowerFitdBm = mag2db(abs(yout)) + 30 - 20;

%去除噪点
inputPowerRange = 20;
idxToDiscard = paInputPowerdBm < (max(paInputPowerdBm)-inputPowerRange);%去除掉与最大输入功率相差20的点

paInputPowerdBm(idxToDiscard) = [];
paOutputPowerdBm(idxToDiscard) = [];
paOutputPowerFitdBm(idxToDiscard) = [];

plot(paInputPowerdBm, paOutputPowerdBm, 'o', paInputPowerdBm,paOutputPowerFitdBm, '.')
grid on
xlabel('Input Power (dBm)')
ylabel('Output Power (dBm)')
title('AM/AM')

%AM/PM图绘制
paInputPhase = angle(x);
paOutputPhase = angle(y);
paOutputPhaseFit = angle(yout);

paInputPhase(idxToDiscard) = [];
paOutputPhase(idxToDiscard) = [];
paOutputPhaseFit(idxToDiscard) = [];

paPhaseChange =  paInputPhase - paOutputPhase;
paPhaseChangeFit =paOutputPhaseFit - paOutputPhase;

%将角度集中在-pi-pi之间
lambdaWrapped = wrapToPi(paPhaseChange);
lambdaWrapped2 = wrapToPi(paPhaseChangeFit);

plot(paInputPowerdBm, lambdaWrapped, 'o', paInputPowerdBm,lambdaWrapped2, '.')
grid on
xlabel('Input Power (dBm)')
ylabel('Phase Change')
title('AM/PM')
```



## 我的程序

### 1、毕设论文1_矩阵处理函数

#### 1、改进后模型

```matlab
function [Xout] = MatrixDeal(x,M,K)


%构造x(n - m)矩阵，m从1到M, 矩阵大小为xLength*M,有效长度为Xlength-M
xnm = x;
xhalfLength = length(xnm);
for m = 1 : M
    xnm(M+1:xhalfLength , m) = x(M-m+1:xhalfLength-m , 1) ;  %(m从1-M)
end

xnm_fix = xnm; %保存x(n - m)矩阵
xnm_fixabs = abs(xnm_fix);%保存 |x(n - m)| 矩阵

xn_fix = x;%保存 x(n) 矩阵
xn_fixabs = abs(xn_fix);%保存 |x(n)| 矩阵

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%m从0开始到M
xnm_0 = x;
xhalfLength = length(xnm_0);
for m = 0 : M
    xnm_0(M+1:xhalfLength , m+1) = x(M-m+1:xhalfLength-m , 1) ;  %(m从1-M)
end
xnm_fix_0 = xnm_0; %保存x(n - m)矩阵
xnm_fixabs_0 = abs(xnm_fix_0);%保存 |x(n - m)| 矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%构造Bk矩阵
Bk = (0:K-1);
Bk = Bk/K;


%f0
%构造 x(n-m)矩阵，矩阵大小为xLength * (M*k)
% X1 = xnm_fix(M+1:xhalfLength,:);

X1 = xnm_fix_0(M+1:xhalfLength,:);

%F21
%构造| ( |x(n-m)| - Bk ) | * x(n - m) * |x(n)| 矩阵，矩阵大小为xLength * (M*k)
xmk_F21 = abs( xnm_fixabs - Bk(1) ) .* xnm_fix .* xn_fixabs;
for k = 2 : K
    mid = abs( xnm_fixabs - Bk(k) ) .* xnm_fix .* xn_fixabs;
    xmk_F21 = cat(2,xmk_F21,mid);
end
%矩阵大小为(xLength-M) * (M*K)
X2 = xmk_F21(M+1:xhalfLength,:);


%f22
%构造 | ( |x(n-m)| - Bk ) | * x(n) 矩阵，矩阵大小为xLength * (M*k)
xmk_f22 = abs( xnm_fixabs - Bk(1) ) .* xn_fix;
for k = 2 : K
    mid = abs( xnm_fixabs - Bk(k) ) .* xn_fix;
    xmk_f22 = cat(2,xmk_f22,mid);
end
%矩阵大小为(xLength-M) * (M*K)
X3 = xmk_f22(M+1:xhalfLength,:);


%F23
%构造  | ( |x(n-m)| - Bk ) | * x(n-m) 矩阵，矩阵大小为xLength * (M*k)
% xmk_F23 = abs( xnm_fixabs - Bk(1) ) .* xnm_fix;
% for k = 2 : K
%     mid = abs( xnm_fixabs - Bk(k) ) .* xnm_fix;
%     xmk_F23 = cat(2,xmk_F23,mid);
% end
% %矩阵大小为(xLength-M) * (M*K)
% X4 = xmk_F23(M+1:xhalfLength,:);



xmk_F23 = abs( xnm_fixabs_0 - Bk(1) ) .* xnm_fix_0;
for k = 2 : K
    mid = abs( xnm_fixabs_0 - Bk(k) ) .* xnm_fix_0;
    xmk_F23 = cat(2,xmk_F23,mid);
end
%矩阵大小为(xLength-M) * (M*K)
X4 = xmk_F23(M+1:xhalfLength,:);


%f24
%构造 | ( |x(n)| - Bk ) | * x(n-m) 矩阵，矩阵大小为xLength * (M*k)
xmk_f24 = abs( xn_fixabs - Bk(1) ) .* xnm_fix;
for k = 2 : K
    mid = abs( xn_fixabs - Bk(k) ) .* xnm_fix;
    xmk_f24 = cat(2,xmk_f24,mid);
end
%矩阵大小为(xLength-M) * (M*K)
X5 = xmk_f24(M+1:xhalfLength,:);


%合并矩阵，矩阵大小为(xLength-M) *（M*K*4+M） 
Xout = [X1 X2 X3 X4 X5];

end
```

#### 2、DVR模型

```matlab
function [Xout] = MatrixDeal_DVR(x,M,K)


%构造x(n - m)矩阵，m从1到M, 矩阵大小为xLength*M,有效长度为Xlength-M-1
xnm = x;
xhalfLength = length(xnm);
for m = 1 : M
    xnm(M+1:xhalfLength , m) = x(M-m+1:xhalfLength-m , 1) ;
end

xnm_fix = xnm; %保存x(n - m)矩阵
xnm_fixabs = abs(xnm_fix);%保存 |x(n - m)| 矩阵
    
xn_fix = x;%保存 x(n) 矩阵
xn_fixabs = abs(xn_fix);%保存 |x(n)| 矩阵

%构造Bk矩阵
Bk = (0:K-1);
Bk = Bk/K;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%m从0开始到M
xnm_0 = x;
xhalfLength = length(xnm_0);
for m = 0 : M
    xnm_0(M+1:xhalfLength , m+1) = x(M-m+1:xhalfLength-m , 1) ;  %(m从1-M)
end
xnm_fix_0 = xnm_0; %保存x(n - m)矩阵
xnm_fixabs_0 = abs(xnm_fix_0);%保存 |x(n - m)| 矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %构造xAngle(n-m)矩阵
% xAngle = angle(xnm);
% for m = 1 : M
%     xAngle(M+1:xhalfLength , m) = xAngle(M-m+1:xhalfLength-m , 1) ;
% end
% xAngle_fix = 1i * xAngle;%保存i * xAngle(n-m)矩阵

%构造xAngle(n-m)矩阵
xAngle = angle(xnm);
for m = 0 : M
    xAngle(M+1:xhalfLength , m+1) = xAngle(M-m+1:xhalfLength-m , 1) ;
end
xAngle_fix = 1i * xAngle;%保存i * xAngle(n-m)矩阵

%f0
%构造 x(n-m)矩阵，矩阵大小为xLength * (M*k)
% X1 = xnm_fix(M+1:xhalfLength,:);

X1 = xnm_fix_0(M+1:xhalfLength,:);

% %f1
% %构造 | ( |x(n-m)| - Bk ) | * e^(xAngle_fix) 矩阵
% xmk_f1 = abs( xnm_fixabs - Bk(1) ) .* exp(xAngle_fix);
% for k = 2 : K
%     mid = abs( xnm_fixabs - Bk(k) ) .* exp(xAngle_fix);
%     xmk_f1 = cat(2,xmk_f1,mid);
% end
% X2 = xmk_f1(M+1:xhalfLength,:);

%f1
%构造 | ( |x(n-m)| - Bk ) | * e^(xAngle_fix) 矩阵
xmk_f1 = abs( xnm_fixabs_0 - Bk(1) ) .* exp(xAngle_fix);
for k = 2 : K
    mid = abs( xnm_fixabs_0 - Bk(k) ) .* exp(xAngle_fix);
    xmk_f1 = cat(2,xmk_f1,mid);
end
X2 = xmk_f1(M+1:xhalfLength,:);


% %f21
% %构造| ( |x(n-m)| - Bk ) | * e^(xAngle_fix) * |x(n)| 矩阵，矩阵大小为xLength * (M*k)
% xmk_f21 = abs( xnm_fixabs - Bk(1) ) .* exp(xAngle_fix) .* xn_fixabs;
% for k = 2 : K
%     mid = abs( xnm_fixabs - Bk(k) ) .* exp(xAngle_fix) .* xn_fixabs;
%     xmk_f21 = cat(2,xmk_f21,mid);
% end
% %矩阵大小为(xLength-M) * (M*K)
% X3 = xmk_f21(M+1:xhalfLength,:);

%f21
%构造| ( |x(n-m)| - Bk ) | * e^(xAngle_fix) * |x(n)| 矩阵，矩阵大小为xLength * (M*k)
xmk_f21 = abs( xnm_fixabs_0 - Bk(1) ) .* exp(xAngle_fix) .* xn_fixabs;
for k = 2 : K
    mid = abs( xnm_fixabs_0 - Bk(k) ) .* exp(xAngle_fix) .* xn_fixabs;
    xmk_f21 = cat(2,xmk_f21,mid);
end
%矩阵大小为(xLength-M) * (M*K)
X3 = xmk_f21(M+1:xhalfLength,:);


%f22
%构造| ( |x(n-m)| - Bk ) | * x(n) 矩阵，矩阵大小为xLength * (M*k)
xmk_f22 = abs( xnm_fixabs - Bk(1) ) .* xn_fix;
for k = 2 : K
    mid = abs( xnm_fixabs - Bk(k) ) .* xn_fix;
    xmk_f22 = cat(2,xmk_f22,mid);
end
%矩阵大小为(xLength-M) * (M*K)
X4 = xmk_f22(M+1:xhalfLength,:);


%f23
%构造| ( |x(n-m)| - Bk ) | * x(n-m) 矩阵，矩阵大小为xLength * (M*k)
xmk_f23 = abs( xnm_fixabs - Bk(1) ) .* xnm_fix;
for k = 2 : K
    mid = abs( xnm_fixabs - Bk(k) ) .* xnm_fix;
    xmk_f23 = cat(2,xmk_f23,mid);
end
%矩阵大小为(xLength-M) * (M*K)
X5 = xmk_f23(M+1:xhalfLength,:);


%f24
%构造 | ( |x(n-m)| - Bk ) | * x(n-m) 矩阵，矩阵大小为xLength * (M*k)
xmk_f24 = abs( xn_fixabs - Bk(1) ) .* xnm_fix;
for k = 2 : K
    mid = abs( xn_fixabs - Bk(k) ) .* xnm_fix;
    xmk_f24 = cat(2,xmk_f24,mid);
end
%矩阵大小为(xLength-M) * (M*K)
X6 = xmk_f24(M+1:xhalfLength,:);

%合并矩阵，矩阵大小为(xLength-M) *（M*K*4+M） 
Xout = [X1 X2 X3 X4 X5 X6];


end
```

#### 3、记忆多项式

```matlab
function [Xout] = MatrixDeal_MP(x,M,K)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明


%构造(x(n) - m)矩阵，m从0到M-1, 矩阵大小为xLength*M
xm = x;
xmLength = length(xm);
for m = 0 : M - 1
    xm(M:xmLength , m+1) = xm(M-m:xmLength-m , 1) ;%有效数据从第M行开始
end
xm_fix = xm; %保存(x(n) - m)矩阵

%构造(x(n) - m) * |x(n) - m|^k 矩阵，矩阵大小为(xLength-M+1)*（M*K）
xmAbs = abs(xm_fix);%保存|x - m|矩阵
for k = 1 : K-1
    mid = (xmAbs.^ k) .* xm_fix;
    xm = cat(2, xm, mid);
end
xmk = xm(M:xmLength,:);

%矩阵大小为(xLength-M) * (M*K)
Xout = xmk;


end
```



### 2、毕设论文1_绘图函数

```matlab
function   paCharPlot(paInput,paOutput,paOutputFit,Type)

paInputMagnitude = abs(paInput) ;
paOutputMagnitude = abs(paOutput) ;
paOutputMagnitudeFit = abs(paOutputFit) ;

idxToDiscard = paOutputMagnitudeFit > 1 | paOutputMagnitude > 1;
paInputMagnitude(idxToDiscard) = [];
paOutputMagnitude(idxToDiscard) = [];
paOutputMagnitudeFit(idxToDiscard) = [];

paGain = paOutputMagnitude - paInputMagnitude;
paGainFit = paOutputMagnitudeFit - paInputMagnitude;


%AM/PM图绘制
paInputPhase = angle(paInput);
paOutputPhase = angle(paOutput);
paOutputPhaseFit = angle(paOutputFit);

paInputPhase(idxToDiscard) = [];
paOutputPhase(idxToDiscard) = [];
paOutputPhaseFit(idxToDiscard) = [];

paPhaseChange =  paInputPhase - paOutputPhase;
paPhaseChangeFit =paInputPhase - paOutputPhaseFit;

%将角度数据集中在-pi-pi之间
lambdaWrapped = wrapToPi(paPhaseChange);
lambdaWrappedFit = wrapToPi(paPhaseChangeFit);    

% lambdaWrapped = rad2deg(lambdaWrapped);
% lambdaWrapped2 = rad2deg(lambdaWrapped2);


switch Type
    case 'AM/AM'
        figure;
        plot(paInputMagnitude, paOutputMagnitude, '.', paInputMagnitude,paOutputMagnitudeFit, '.')
        grid on
        xlabel('Input Magnitude')
        ylabel('Output Magnitude')
        legend({'Actual','Model/Linear'},'Location','northwest')
        title('AM/AM')


    case 'PM/AM'
        figure;
        plot(paInputMagnitude, lambdaWrapped, '.', paInputMagnitude,lambdaWrappedFit, '.')
        grid on
        xlabel('Input Magnitude')
        ylabel('Phase Change')
        legend({'Actual','Model/Linear'},'Location','northwest')
        title('AM/PM')

       
    case 'Gain'
        figure;
        plot(paInputMagnitude, paGain, '.', paInputMagnitude,paGainFit, '.')
        grid on
        xlabel('Input Magnitude')
        ylabel('Gain')
        legend({'Actual','Model/Linear'},'Location','northwest')
        title('PA Gain')


    case 'Couple'
        yyaxis left;
        plot(paInputPowerdBm, paOutputPowerdBm, 'o', paInputPowerdBm,paOutputPowerFitdBm, '.');
        yyaxis right;
        plot(paInputPowerdBm, lambdaWrapped, 'o', paInputPowerdBm,lambdaWrapped2, '.');
        hold on

end


end
```

### 3、毕设论文1_模型建立

#### 1、改进后模型

```matlab
%本例使用NXP Airfast LDMOS Doherty PA，工作频率3.6-3.8 GHz，增益29 dB。
%带宽100MHZ    100000000
%采样率sampleRate = 860160000
%信号类型OFDM
%If testSignal is "OFDM", this example uses a 5G-like OFDM waveform with 64-QAM modulated signals for each subcarrier.
% If testSignal is "Tones", this example uses two tones at 1.8 MHz and 2.6 MHz, to test the intermodulation caused by the PA.


sampleRate = 860160000;
% sampleRate = 430080000;
testSignal = 'OFDM';

%读取数据
% x = (0 : 5);
x = paInput;
% x = txData1;
% x = paInput_40MHZ;
x = x(:);
xLength = length(x);
%使用前半部分数据用于参数提取，后半部分数据用于模型验证
half = round(xLength/2);

% y = (2 : 7);
y = paOutput;
% y = tout;
% y = paOutput_40MHZ;
y = y(:);
yLength = length(y);

%将数据归一化
% x = x / abs(max(x));
% y = y / abs(max(y));

M = 2;
K = 8;

%计算PA系数矩阵，矩阵大小为（M*K*4+M） * 1

%使用前一半数据估计参数，将数据归一化
XcoefPA = x(1:half) / abs ( max(x(1:half)) );  
YcoefPA = y(1:half) / abs ( max(y(1:half)) );  

% XcoefPA = x(1:half);  
% YcoefPA = y(1:half);  

%去除多余数据
YcoefPA = YcoefPA(M+1:half) ;                   %经过处理后，数据去除掉前M个
coefPA = MatrixDeal(XcoefPA,M,K) \ YcoefPA;

%最小二乘算法
% coefPA  = inv( (XcoefPA') *XcoefPA ) * (XcoefPA') * YcoefPA;
% coefPA  = ( (XcoefPA') *XcoefPA ) \ (XcoefPA') * YcoefPA;

%带入验证PA模型输出，矩阵大小为xLength*1
XmodelPA = x(half:xLength) / abs ( max(x(half:xLength)) );      %使用后一半的数据验证模型，将数据归一化
YmodelPA = y(half:xLength) / abs ( max(y(half:xLength)) ); 

% XmodelPA = x(half:xLength);      %使用后一半的数据验证模型，将数据归一化
% YmodelPA = y(half:xLength); 

YmodelFitPA = MatrixDeal(XmodelPA,M,K) * coefPA;   %计算模型输出（输出后的数据减少前M个）

%去除多余数据
XmodelPA = XmodelPA(M+1:length(XmodelPA));
YmodelPA = YmodelPA(M+1:length(YmodelPA));

%通过归一化均方误差衡量功放的建模精度
NMSE = 10 * log10( sum ( (abs(YmodelPA - YmodelFitPA)).^2 ) / sum ( (abs(YmodelPA)).^2 ) ) ;

%计算未线性化的EVM
EVM_withoutDPD = sqrt (  sum ( (abs(YmodelPA - XmodelPA)).^2 ) / sum ( (abs(XmodelPA)).^2) ) * 100 ;
disp(['The EVM_withoutDPD is ' num2str(EVM_withoutDPD) '%']);

%绘图
paCharPlot(XmodelPA,YmodelPA,YmodelFitPA,'AM/AM');
% paCharPlot(XmodelPA,YmodelPA,XmodelPA,'AM/AM');
paCharPlot(XmodelPA,YmodelPA,YmodelFitPA,'PM/AM');
% paCharPlot(XmodelPA,YmodelPA,YmodelFitPA,'Gain');


%观察频谱
%sa = helperPACharPlotSpectrum([YmodelPA  YmodelFitPA],{'Actual PA Output', 'Model PA Output'}, sampleRate,testSignal);

%预失真器建立

%预失真器参数提取
XcoefDPD = x(1:half) / abs ( max(x(1:half)) );  %使用前一半数据估计参数，将数据归一化
YcoefDPD = y(1:half) / abs ( max(y(1:half)) );  

% XcoefDPD = x(1:half);  %使用前一半数据估计参数，将数据归一化
% YcoefDPD = y(1:half);  

XcoefDPD = XcoefDPD(M+1:half);                     %将PA实际输入数据作为输出数据，取前一半数据
coefDPD = MatrixDeal(YcoefDPD,M,K) \ XcoefDPD;      %计算得到DPD参数模型，与PA行为模型一致

%计算DPD模型输出
XmodelDPD = x(half:xLength) / abs ( max(x(half:xLength)) );         %使用后一半数据进行DPD模型验证
YmodelDPD = y(half:xLength) / abs ( max(y(half:xLength)) ); 

% XmodelDPD = x(half:xLength);         %使用后一半数据进行DPD模型验证
% YmodelDPD = y(half:xLength); 

YmodelFitDPD = MatrixDeal(XmodelDPD,M,K) * coefDPD;  %计算得到DPD输出（输出数据减少M个）


%去除多余数据
% %对比PA模型与DPD模型的非线性特性
% XmodelDPD = x(half+M:xLength,:);         %实际用于DPD模型验证的输入数据(去除前M个)
% YmodelDPD = y(half+M:yLength,:);         %实际的DPD模型的输出数据（去除前M个）
% paCharPlot(XmodelDPD,YmodelFitPA,YmodelFitDPD,'AM/AM');


%将DPD模型输出，输入到PA模型中
Ylinear = MatrixDeal(YmodelFitDPD,M,K) * coefPA;

%计算误差向量幅度EVM
XDPDPA = XmodelDPD(M+M+1:length(XmodelDPD));
YDPDPA = YmodelDPD(M+M+1:length(YmodelDPD));

EVM_withDPD = sqrt (  sum ( (abs(Ylinear - XDPDPA)).^2 ) / sum ( (abs(XDPDPA)).^2) ) * 100 ;
disp(['The EVM_withDPD is ' num2str(EVM_withDPD) '%']);

%绘图
paCharPlot(XDPDPA,YDPDPA,Ylinear,'AM/AM');
% % paCharPlot(XDPDPA,XDPDPA,Ylinear,'AM/AM');
paCharPlot(XDPDPA,YDPDPA,Ylinear,'PM/AM');
% paCharPlot(XDPDPA,YDPDPA,Ylinear,'Gain');

%邻信道功率比ACPR
%sa = helperPACharPlotSpectrum([YDPDPA  Ylinear],{'Actual PA Output', 'with DPD Output'}, sampleRate,testSignal);


```

#### 2、DVR模型

```matlab
%本例使用NXP Airfast LDMOS Doherty PA，工作频率3.6-3.8 GHz，增益29 dB。
%带宽100MHZ
%采样率sampleRate = 860160000
%信号类型OFDM
%If testSignal is "OFDM", this example uses a 5G-like OFDM waveform with 64-QAM modulated signals for each subcarrier.
% If testSignal is "Tones", this example uses two tones at 1.8 MHz and 2.6 MHz, to test the intermodulation caused by the PA.


sampleRate = 860160000;
% sampleRate = 430080000;
testSignal = 'OFDM';

%读取数据
% x = (0 : 5);
x = paInput;
% x = txData1;
% x = paInput_40MHZ;
x = x(:);
xLength = length(x);
%使用前半部分数据用于参数提取，后半部分数据用于模型验证
half = round(xLength/2);

% y = (2 : 7);
y = paOutput;
% y = tout;
% y = paOutput_40MHZ;
y = y(:);
yLength = length(y);

%将数据归一化
% x = x / abs(max(x));
% y = y / abs(max(y));



M = 2;
K = 8;

%计算PA系数矩阵，矩阵大小为（M*K*4+M） * 1

%使用前一半数据估计参数，将数据归一化
XcoefPA = x(1:half) / abs ( max(x(1:half)) );  
YcoefPA = y(1:half) / abs ( max(y(1:half)) );  

% XcoefPA = x(1:half);  
% YcoefPA = y(1:half);  

%去除多余数据
YcoefPA = YcoefPA(M+1:half) ;                   %经过处理后，数据去除掉前M个
coefPA = MatrixDeal_DVR(XcoefPA,M,K) \ YcoefPA;

%最小二乘算法
% coefPA  = inv( (XcoefPA') *XcoefPA ) * (XcoefPA') * YcoefPA;
% coefPA  = ( (XcoefPA') *XcoefPA ) \ (XcoefPA') * YcoefPA;

%带入验证PA模型输出，矩阵大小为xLength*1
XmodelPA = x(half:xLength) / abs ( max(x(half:xLength)) );      %使用后一半的数据验证模型，将数据归一化
YmodelPA = y(half:xLength) / abs ( max(y(half:xLength)) ); 

% XmodelPA = x(half:xLength);      %使用后一半的数据验证模型，将数据归一化
% YmodelPA = y(half:xLength); 

YmodelFitPA = MatrixDeal_DVR(XmodelPA,M,K) * coefPA;   %计算模型输出（输出后的数据减少前M个）

%去除多余数据
XmodelPA = XmodelPA(M+1:length(XmodelPA));
YmodelPA = YmodelPA(M+1:length(YmodelPA));

%通过归一化均方误差衡量功放的建模精度
NMSE = 10 * log10( sum ( (abs(YmodelPA - YmodelFitPA)).^2 ) / sum ( (abs(YmodelPA)).^2 ) ) ;

%计算未线性化的EVM
EVM_withoutDPD = sqrt (  sum ( (abs(YmodelPA - XmodelPA)).^2 ) / sum ( (abs(XmodelPA)).^2) ) * 100 ;
disp(['The EVM_withoutDPD is ' num2str(EVM_withoutDPD) '%']);

%绘图
paCharPlot(XmodelPA,YmodelPA,YmodelFitPA,'AM/AM');
% paCharPlot(XmodelPA,YmodelPA,XmodelPA,'AM/AM');
paCharPlot(XmodelPA,YmodelPA,YmodelFitPA,'PM/AM');


%观察频谱
%sa = helperPACharPlotSpectrum([YmodelPA  YmodelFitPA],{'Actual PA Output', 'Model PA Output'}, sampleRate,testSignal);


%预失真器建立

%预失真器参数提取
XcoefDPD = x(1:half) / abs ( max(x(1:half)) );  %使用前一半数据估计参数，将数据归一化
YcoefDPD = y(1:half) / abs ( max(y(1:half)) );  

% XcoefDPD = x(1:half);  %使用前一半数据估计参数，将数据归一化
% YcoefDPD = y(1:half);  

XcoefDPD = XcoefDPD(M+1:half);                     %将PA实际输入数据作为输出数据，取前一半数据
coefDPD = MatrixDeal_DVR(YcoefDPD,M,K) \ XcoefDPD;      %计算得到DPD参数模型，与PA行为模型一致

%计算DPD模型输出
XmodelDPD = x(half:xLength) / abs ( max(x(half:xLength)) );         %使用后一半数据进行DPD模型验证
YmodelDPD = y(half:xLength) / abs ( max(y(half:xLength)) ); 

% XmodelDPD = x(half:xLength);         %使用后一半数据进行DPD模型验证
% YmodelDPD = y(half:xLength); 

YmodelFitDPD = MatrixDeal_DVR(XmodelDPD,M,K) * coefDPD;  %计算得到DPD输出（输出数据减少M个）


%去除多余数据
% %对比PA模型与DPD模型的非线性特性
% XmodelDPD = x(half+M:xLength,:);         %实际用于DPD模型验证的输入数据(去除前M个)
% YmodelDPD = y(half+M:yLength,:);         %实际的DPD模型的输出数据（去除前M个）
% paCharPlot(XmodelDPD,YmodelFitPA,YmodelFitDPD,'AM/AM');


%将DPD模型输出，输入到PA模型中
Ylinear = MatrixDeal_DVR(YmodelFitDPD,M,K) * coefPA;

%计算误差向量幅度EVM
XDPDPA = XmodelDPD(M+M+1:length(XmodelDPD));
YDPDPA = YmodelDPD(M+M+1:length(YmodelDPD));

EVM_withDPD = sqrt (  sum ( (abs(Ylinear - XDPDPA)).^2 ) / sum ( (abs(XDPDPA)).^2) ) * 100 ;
disp(['The EVM_withDPD is ' num2str(EVM_withDPD) '%']);

%绘图
paCharPlot(XDPDPA,YDPDPA,Ylinear,'AM/AM');
% paCharPlot(XDPDPA,XDPDPA,Ylinear,'AM/AM');
paCharPlot(XDPDPA,YDPDPA,Ylinear,'PM/AM');

%邻信道功率比ACPR
%sa = helperPACharPlotSpectrum([YDPDPA  Ylinear],{'Actual PA Output', 'with DPD Output'}, sampleRate,testSignal);


```

#### 3、MP模型

```matlab
% % 为了表征AM/AM传递函数，计算输入功率值范围内的平均输出功率。
% % 测量的单位是伏，总体阻抗为100欧姆，
% % 在发射器和接收器之间划分。
% % 将测量的基带样本转换为dBm的功率值。
% % +30 dB项用于dBW到dBm的转换，
% % -20 dB项用于100欧姆阻抗。
% 
% %mag2db函数用于将振幅转换为dB 
% %原始输入输出数据是复数，绝对值转换为功率
% paInputdBm  = mag2db(abs(paInput)) + 30 - 20;
% paOutputdBm  = mag2db(abs(paOutput)) + 30 - 20;


% modType = 'memPoly';
% 
% memLen = 3;
% degLen = 5;
% 
% numDataPts = length(paInput);
% halfDataPts = round(numDataPts/2);
% 
% fitCoefMatMem = helperPACharMemPolyModel('coefficientFinder', paInput(1:halfDataPts),paOutput(1:halfDataPts),memLen,degLen,modType);
% 
% disp(abs(fitCoefMatMem));
% 
% rmsErrorTimeMem = helperPACharMemPolyModel('errorMeasure', paInput, paOutput, fitCoefMatMem, modType);
% 
% disp(['Percent RMS error in time domain is ' num2str(rmsErrorTimeMem) '%']);
% paOutputFitMem = helperPACharMemPolyModel('signalGenerator',  paInput, fitCoefMatMem, modType);
%   
% helperPACharPlotTime(paOutput, paOutputFitMem, sampleRate);
% helperPACharPlotGain(paInput, paOutput, paOutputFitMem);




%本例使用NXP Airfast LDMOS Doherty PA，工作频率3.6-3.8 GHz，增益29 dB。
%带宽100MHZ
%采样率sampleRate = 860160000
%信号类型OFDM
%If testSignal is "OFDM", this example uses a 5G-like OFDM waveform with 64-QAM modulated signals for each subcarrier.
% If testSignal is "Tones", this example uses two tones at 1.8 MHz and 2.6 MHz, to test the intermodulation caused by the PA.


sampleRate = 860160000;
% sampleRate = 430080000;
testSignal = 'OFDM';

%读取数据
% x = (0 : 5);
x = paInput;
% x = txData1;
% x = paInput_40MHZ;
x = x(:);
xLength = length(x);
%使用前半部分数据用于参数提取，后半部分数据用于模型验证
half = round(xLength/2);

% y = (2 : 7);
y = paOutput;
% y = tout;
% y = paOutput_40MHZ;
y = y(:);
yLength = length(y);

%将数据归一化
% x = x / abs(max(x));
% y = y / abs(max(y));

M = 2;
K = 8;


%使用前一半数据用于PA参数提取
XcoefPA = x(1:half) / abs ( max(x(1:half)) );  
YcoefPA = y(1:half) / abs ( max(y(1:half)) );  

% XcoefPA = x(1:half);  
% YcoefPA = y(1:half);  

YcoefPA = YcoefPA(M:half) ;                   %经过处理后，数据去除掉前M个
coefPA = MatrixDeal_MP(XcoefPA,M,K) \ YcoefPA;

%使用后一半数据用于PA模型验证
XmodelPA = x(half:xLength) / abs ( max(x(half:xLength)) );      %使用后一半的数据验证模型，将数据归一化
YmodelPA = y(half:xLength) / abs ( max(y(half:xLength)) ); 
% XmodelPA = x(half:xLength);      %使用后一半的数据验证模型，将数据归一化
% YmodelPA = y(half:xLength); 

YmodelFitPA = MatrixDeal_MP(XmodelPA,M,K) * coefPA;   %计算模型输出（输出后的数据减少前M个）

%去除多余数据
XmodelPA = XmodelPA(M:length(XmodelPA));
YmodelPA = YmodelPA(M:length(YmodelPA));

%计算未线性化的EVM
EVM_withoutDPD = sqrt (  sum ( (abs(YmodelPA - XmodelPA)).^2 ) / sum ( (abs(XmodelPA)).^2) ) * 100 ;
disp(['The EVM_withoutDPD is ' num2str(EVM_withoutDPD) '%']);

%计算NMSE
NMSE = 10 * log10( sum ( (abs(YmodelPA - YmodelFitPA)).^2 ) / sum ( (abs(YmodelPA)).^2 ) ) ;

%绘图
paCharPlot(XmodelPA,YmodelPA,YmodelFitPA,'AM/AM');
paCharPlot(XmodelPA,YmodelPA,YmodelFitPA,'PM/AM');



%使用前一半数据估计参数，将数据归一化
XcoefDPD = x(1:half) / abs ( max(x(1:half)) );  %使用前一半数据估计参数，将数据归一化
YcoefDPD = y(1:half) / abs ( max(y(1:half)) );  

% XcoefDPD = x(1:half);  %使用前一半数据估计参数，将数据归一化
% YcoefDPD = y(1:half);  


XcoefDPD = XcoefDPD(M:half);                     %将PA实际输入数据作为输出数据，取前一半数据
coefDPD = MatrixDeal_MP(YcoefDPD,M,K) \ XcoefDPD;      %计算得到DPD参数模型，与PA行为模型一致


%模型验证
XmodelDPD = x(half:xLength) / abs ( max(x(half:xLength)) );         %使用后一半数据进行DPD模型验证
YmodelDPD = y(half:xLength) / abs ( max(y(half:xLength)) ); 

% XmodelDPD = x(half:xLength);         %使用后一半数据进行DPD模型验证
% YmodelDPD = y(half:xLength); 

YmodelFitDPD = MatrixDeal_MP(XmodelDPD,M,K) * coefDPD;  %计算得到DPD输出（输出数据减少M个）

Ylinear = MatrixDeal_MP(YmodelFitDPD,M,K) * coefPA;

%计算误差向量幅度EVM
XDPDPA = XmodelDPD(M+M-1:length(XmodelDPD));
YDPDPA = YmodelDPD(M+M-1:length(YmodelDPD));

EVM_withDPD = sqrt (  sum ( (abs(Ylinear - XDPDPA)).^2 ) / sum ( (abs(XDPDPA)).^2) ) * 100 ;
disp(['The EVM_withDPD is ' num2str(EVM_withDPD) '%']);

%绘图
paCharPlot(XDPDPA,YDPDPA,Ylinear,'AM/AM');
paCharPlot(XDPDPA,YDPDPA,Ylinear,'PM/AM');

%邻信道功率比ACPR
sa = helperPACharPlotSpectrum([YDPDPA  Ylinear],{'Actual PA Output', 'with DPD Output'}, sampleRate,testSignal);

```

