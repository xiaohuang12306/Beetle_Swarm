%% 清空环境
clc;clear all;close all;
%% 目标函数
%    下列函数为：
%    y=0.072*x1+0.063*x2+0.057*x3+0.05*x4+0.032*x5+0.0442*x6+0.0675*x7+7
%    约束条件为：
%    0.072*x1+0.063*x2+0.057*x3+0.05*x4+0.032*x5+0.0442*x6+0.0675*x7<=264.4
%    128*x1+78.1*x2+64.1*x3+43*x4+58.1*x5+36.9*x6+50.5*x7<=67919
a=[0.072,0.063,0.057,0.05,0.032,0.0442,0.0675];
b=[0.072,0.063,0.057,0.005,0.032,0.0442,0.0675];
c=[128,78.1,64.1,43,58.1,36.9,50.5];
fun= @(X)(a*X+7);
cons1= @(X)(b*X<=264.4);
cons2= @(X)(c*X>=-20);
%% 设置种群参数
%   需要自行配置
sizepop = 500;                         % 初始种群个数
dim = 7;                           % 空间维数
ger = 500;                       % 最大迭代次数    
xlimit_max = 50*ones(dim,1);     % 设置位置参数限制(矩阵的形式可以多维)
xlimit_min = -50*ones(dim,1);
vlimit_max = 1*ones(dim,1);      % 设置速度限制
vlimit_min = -1*ones(dim,1);
c_1 = 0.8;                       % 惯性权重
c_2 = 0.5;                       % 自我学习因子
c_3 = 0.5;                       % 群体学习因子 

%% 生成初始种群
%  首先随机生成初始种群位置
%  然后随机生成初始种群速度
%  然后初始化个体历史最佳位置，以及个体历史最佳适应度
%  然后初始化群体历史最佳位置，以及群体历史最佳适应度
for i=1:dim
    for j=1:sizepop
        pop_x(i,j) = xlimit_min(i)+(xlimit_max(i) - xlimit_min(i))*rand;  % 初始种群的位置
        pop_v(i,j) = vlimit_min(i)+(vlimit_max(i) - vlimit_min(i))*rand;  % 初始种群的速度
    end
end                 
gbest = pop_x;                                % 每个个体的历史最佳位置
for j=1:sizepop
    if cons1(pop_x(:,j))
        if cons2(pop_x(:,j))
            fitness_gbest(j) = fun(pop_x(:,j));                      % 每个个体的历史最佳适应度
        else
            fitness_gbest(j) = 10^10; 
        end
    else
        fitness_gbest(j) = 10^10; 
    end  
end                  
zbest = pop_x(:,1);                           % 种群的历史最佳位置
fitness_zbest = fitness_gbest(1);             % 种群的历史最佳适应度
for j=1:sizepop
    if fitness_gbest(j) < fitness_zbest       % 如果求最小值，则为<; 如果求最大值，则为>; 
        zbest = pop_x(:,j);
        fitness_zbest=fitness_gbest(j);
    end
end


%% 粒子群迭代
%    更新速度并对速度进行边界处理    
%    更新位置并对位置进行边界处理
%    进行自适应变异
%    进行约束条件判断并计算新种群各个个体位置的适应度
%    新适应度与个体历史最佳适应度做比较
%    个体历史最佳适应度与种群历史最佳适应度做比较
%    再次循环或结束

iter = 1;                        %迭代次数
record = zeros(ger, 1);          % 记录器
while iter <= ger
    for j=1:sizepop
        %    更新速度并对速度进行边界处理 
        pop_v(:,j)= c_1 * pop_v(:,j) + c_2*rand*(gbest(:,j)-pop_x(:,j))+c_3*rand*(zbest-pop_x(:,j));% 速度更新
        for i=1:dim
            if  pop_v(i,j) > vlimit_max(i)
                pop_v(i,j) = vlimit_max(i);
            end
            if  pop_v(i,j) < vlimit_min(i)
                pop_v(i,j) = vlimit_min(i);
            end
        end
        
        %    更新位置并对位置进行边界处理
        pop_x(:,j) = pop_x(:,j) + pop_v(:,j);% 位置更新
        for i=1:dim
            if  pop_x(i,j) > xlimit_max(i)
                pop_x(i,j) = xlimit_max(i);
            end
            if  pop_x(i,j) < xlimit_min(i)
                pop_x(i,j) = xlimit_min(i);
            end
        end
        
        %    进行自适应变异
        if rand > 0.85
            i=ceil(dim*rand);
            pop_x(i,j)=xlimit_min(i) + (xlimit_max(i) - xlimit_min(i)) * rand;
        end
  
        %    进行约束条件判断并计算新种群各个个体位置的适应度
        if cons1(pop_x(:,j))
            if cons2(pop_x(:,j))
                fitness_pop(j) = fun(pop_x(:,j));                      % 当前个体的适应度
            else
                fitness_pop(j) = 10^10; 
            end
        else
            fitness_pop(j) = 10^10; 
        end
        
        %    新适应度与个体历史最佳适应度做比较
        if fitness_pop(j) < fitness_gbest(j)       % 如果求最小值，则为<; 如果求最大值，则为>; 
            gbest(:,j) = pop_x(:,j);               % 更新个体历史最佳位置            
            fitness_gbest(j) = fitness_pop(j);     % 更新个体历史最佳适应度
        end   
        
        %    个体历史最佳适应度与种群历史最佳适应度做比较
        if fitness_gbest(j) < fitness_zbest        % 如果求最小值，则为<; 如果求最大值，则为>; 
            zbest = gbest(:,j);                    % 更新群体历史最佳位置  
            fitness_zbest=fitness_gbest(j);        % 更新群体历史最佳适应度  
        end    
    end
    
    record(iter) = fitness_zbest;%最大值记录
    
    iter = iter+1;

end
%% 迭代结果输出

plot(record);title('收敛过程')
disp(['最优值：',num2str(fitness_zbest)]);
disp('变量取值：');
disp(zbest);
