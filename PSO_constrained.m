%% ��ջ���
clc;clear all;close all;
%% Ŀ�꺯��
%    ���к���Ϊ��
%    y=0.072*x1+0.063*x2+0.057*x3+0.05*x4+0.032*x5+0.0442*x6+0.0675*x7+7
%    Լ������Ϊ��
%    0.072*x1+0.063*x2+0.057*x3+0.05*x4+0.032*x5+0.0442*x6+0.0675*x7<=264.4
%    128*x1+78.1*x2+64.1*x3+43*x4+58.1*x5+36.9*x6+50.5*x7<=67919
a=[0.072,0.063,0.057,0.05,0.032,0.0442,0.0675];
b=[0.072,0.063,0.057,0.005,0.032,0.0442,0.0675];
c=[128,78.1,64.1,43,58.1,36.9,50.5];
fun= @(X)(a*X+7);
cons1= @(X)(b*X<=264.4);
cons2= @(X)(c*X>=-20);
%% ������Ⱥ����
%   ��Ҫ��������
sizepop = 500;                         % ��ʼ��Ⱥ����
dim = 7;                           % �ռ�ά��
ger = 500;                       % ����������    
xlimit_max = 50*ones(dim,1);     % ����λ�ò�������(�������ʽ���Զ�ά)
xlimit_min = -50*ones(dim,1);
vlimit_max = 1*ones(dim,1);      % �����ٶ�����
vlimit_min = -1*ones(dim,1);
c_1 = 0.8;                       % ����Ȩ��
c_2 = 0.5;                       % ����ѧϰ����
c_3 = 0.5;                       % Ⱥ��ѧϰ���� 

%% ���ɳ�ʼ��Ⱥ
%  ����������ɳ�ʼ��Ⱥλ��
%  Ȼ��������ɳ�ʼ��Ⱥ�ٶ�
%  Ȼ���ʼ��������ʷ���λ�ã��Լ�������ʷ�����Ӧ��
%  Ȼ���ʼ��Ⱥ����ʷ���λ�ã��Լ�Ⱥ����ʷ�����Ӧ��
for i=1:dim
    for j=1:sizepop
        pop_x(i,j) = xlimit_min(i)+(xlimit_max(i) - xlimit_min(i))*rand;  % ��ʼ��Ⱥ��λ��
        pop_v(i,j) = vlimit_min(i)+(vlimit_max(i) - vlimit_min(i))*rand;  % ��ʼ��Ⱥ���ٶ�
    end
end                 
gbest = pop_x;                                % ÿ���������ʷ���λ��
for j=1:sizepop
    if cons1(pop_x(:,j))
        if cons2(pop_x(:,j))
            fitness_gbest(j) = fun(pop_x(:,j));                      % ÿ���������ʷ�����Ӧ��
        else
            fitness_gbest(j) = 10^10; 
        end
    else
        fitness_gbest(j) = 10^10; 
    end  
end                  
zbest = pop_x(:,1);                           % ��Ⱥ����ʷ���λ��
fitness_zbest = fitness_gbest(1);             % ��Ⱥ����ʷ�����Ӧ��
for j=1:sizepop
    if fitness_gbest(j) < fitness_zbest       % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>; 
        zbest = pop_x(:,j);
        fitness_zbest=fitness_gbest(j);
    end
end


%% ����Ⱥ����
%    �����ٶȲ����ٶȽ��б߽紦��    
%    ����λ�ò���λ�ý��б߽紦��
%    ��������Ӧ����
%    ����Լ�������жϲ���������Ⱥ��������λ�õ���Ӧ��
%    ����Ӧ���������ʷ�����Ӧ�����Ƚ�
%    ������ʷ�����Ӧ������Ⱥ��ʷ�����Ӧ�����Ƚ�
%    �ٴ�ѭ�������

iter = 1;                        %��������
record = zeros(ger, 1);          % ��¼��
while iter <= ger
    for j=1:sizepop
        %    �����ٶȲ����ٶȽ��б߽紦�� 
        pop_v(:,j)= c_1 * pop_v(:,j) + c_2*rand*(gbest(:,j)-pop_x(:,j))+c_3*rand*(zbest-pop_x(:,j));% �ٶȸ���
        for i=1:dim
            if  pop_v(i,j) > vlimit_max(i)
                pop_v(i,j) = vlimit_max(i);
            end
            if  pop_v(i,j) < vlimit_min(i)
                pop_v(i,j) = vlimit_min(i);
            end
        end
        
        %    ����λ�ò���λ�ý��б߽紦��
        pop_x(:,j) = pop_x(:,j) + pop_v(:,j);% λ�ø���
        for i=1:dim
            if  pop_x(i,j) > xlimit_max(i)
                pop_x(i,j) = xlimit_max(i);
            end
            if  pop_x(i,j) < xlimit_min(i)
                pop_x(i,j) = xlimit_min(i);
            end
        end
        
        %    ��������Ӧ����
        if rand > 0.85
            i=ceil(dim*rand);
            pop_x(i,j)=xlimit_min(i) + (xlimit_max(i) - xlimit_min(i)) * rand;
        end
  
        %    ����Լ�������жϲ���������Ⱥ��������λ�õ���Ӧ��
        if cons1(pop_x(:,j))
            if cons2(pop_x(:,j))
                fitness_pop(j) = fun(pop_x(:,j));                      % ��ǰ�������Ӧ��
            else
                fitness_pop(j) = 10^10; 
            end
        else
            fitness_pop(j) = 10^10; 
        end
        
        %    ����Ӧ���������ʷ�����Ӧ�����Ƚ�
        if fitness_pop(j) < fitness_gbest(j)       % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>; 
            gbest(:,j) = pop_x(:,j);               % ���¸�����ʷ���λ��            
            fitness_gbest(j) = fitness_pop(j);     % ���¸�����ʷ�����Ӧ��
        end   
        
        %    ������ʷ�����Ӧ������Ⱥ��ʷ�����Ӧ�����Ƚ�
        if fitness_gbest(j) < fitness_zbest        % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>; 
            zbest = gbest(:,j);                    % ����Ⱥ����ʷ���λ��  
            fitness_zbest=fitness_gbest(j);        % ����Ⱥ����ʷ�����Ӧ��  
        end    
    end
    
    record(iter) = fitness_zbest;%���ֵ��¼
    
    iter = iter+1;

end
%% ����������

plot(record);title('��������')
disp(['����ֵ��',num2str(fitness_zbest)]);
disp('����ȡֵ��');
disp(zbest);
