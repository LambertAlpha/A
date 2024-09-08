clc;close all
warning off
luoju=55e-2; % 螺距
k=luoju/2/pi; % 螺线方程的系数 r=k theta
L1=341e-2;
D1=L1-27.5e-2*2; % 龙头把手两个孔之间的距离
L2=220e-2;
D2=L2-27.5e-2*2; % 其他凳子把手两个孔之间的距离

%% 先画出部分螺线
theta=20*2*pi:-0.01:0*pi;
r=k*theta;
x=r.*cos(theta);
y=r.*sin(theta);
figure(1)
set(gcf,'Position',[200 200 600 600]);
plot(x,y,'--')
axis equal
grid on
xlabel('x')
ylabel('y')
hold on

% 第一步，先确定300s内，龙头第一个把手的位置演化
mydtheta=@(t,theta)-1./(k*sqrt(1+theta.^2));
theta0=2*pi*16; % 初始位置时候的角度
dt=0.05; % 时间步长 定为0.05
tspan=0:dt:300; % 求解时间点
[tt,theta]=ode45(mydtheta,tspan,theta0); % 龙格库塔法求解
X1=k*theta.*cos(theta);
Y1=k*theta.*sin(theta);
for i=1:10:length(theta)
    title({['t=',num2str(tt(i))],'龙头前把手的运动轨迹'})
    plot(X1(i),Y1(i),'b.','MarkerSize',10)
    drawnow
end
hwait=waitbar(0,'计算开始...')
%% 第二步,确定每个时间点下，头部凳子的后面一个孔（也要在螺线上）,以及龙身和龙尾凳子各个孔所在位置(都要在螺线上）
N=223; % 龙头+龙身+龙尾总的个数
X=zeros(N+1,length(X1)); 
Y=zeros(N+1,length(X1)); % 每一行代表每个凳子的前把手孔的位置在各个时间点处的值,因为尾部还有一个孔，所以一共223+1个，所以X和Y对应223+1行
Theta=zeros(N+1,length(X1)); % 记录每个孔在各个时刻处的位置对应的角度theta（这为了求速度的）
X(1,:)=X1;
Y(1,:)=Y1;% 第一行已知了，上面求得的头部第一个把手位置数据
Theta(1,:)=theta; % 第一行已知了，第一个把手的角度数据，上面求了
for j=1:length(tt)
    for i=2:N+1 % 在每一个时间点下，对每一行循环计算，意思是求出此时各个把手孔的位置信息
        d=D1*(i<=2)+D2*(i>2); % 分辨下是第一个凳子还是其他凳子，孔之间的距离不一样！
        thetaij=solve_theta(luoju,X(i-1,j),Y(i-1,j),Theta(i-1,j),d); % 子函数求解下一个孔的角度值
        Theta(i,j)=thetaij;
        X(i,j)=k*thetaij*cos(thetaij);
        Y(i,j)=k*thetaij*sin(thetaij);
        waitbar(((j-1)*N+i)/(length(tt)*N),hwait,'已经完成...')
    end
end
close(hwait)
%% 下面来可视化，盘入的动态图
theta=20*2*pi:-0.01:0*pi;
r=k*theta;
x=r.*cos(theta);
y=r.*sin(theta);
set(gcf,'Position',[200 200 600 600])
for j=1:length(tt)
    figure(2)
    plot(x,y,'-')
    axis equal
    grid on
    xlabel('x')
    ylabel('y')
    hold on
    
    plot(X(:,j),Y(:,j),'k-','LineWidth',1.2,'Marker','o','MarkerSize',6,'MarkerFaceColor','b')
    title({['t=',num2str(tt(j))],'盘入的轨迹'})
    drawnow
    hold off
end
%% 下面来求解速度数据(有各个孔的角度数据Theta，速度数据可以用数值差分方法
V=zeros(size(X)); % 每一行对应一个孔的速度在不同时间节点上的值
V(:,1)=-k*sqrt(1+Theta(:,1).^2).*(Theta(:,2)-Theta(:,1))/dt; % 第一个时间点，用前向差分求导数
V(:,end)=-k*sqrt(1+Theta(:,end).^2).*(Theta(:,end)-Theta(:,end-1))/dt; % 最后一个时间点，用后向差分求导数
V(:,2:end-1)=-k*sqrt(1+Theta(:,2:end-1).^2).*(Theta(:,3:end)-Theta(:,1:end-2))/2/dt; % 中间的每一个时间点，用中心差分来求导数
% 不明白为什么这样求导的，去看我的推文：https://mp.weixin.qq.com/s/rI9B-VHT9i2jA28CN7oFMQ 
% 下面验证一下第一个孔位置的速度是不是1，也就是V矩阵的第一行的值
figure
plot(tt,V(1,:),'b-','LineWidth',1.3)
ylim([0 1.1])
xlabel('时间')
ylabel('头把手的速度')
title('验证数值计算得到的速度')

%% %%%%%%%%%%%%%%%%%%%% 到此，我们已经得到了所有把手中心的位置和速度随时间的数据,下面填表
nn=1/dt;
index=1:nn:length(tt); % 找出特定的时间点,来记录
Dataxy=zeros(2*(N+1),length(index)); % 这些记录所有点的位置信息
Dataxy(1:2:end,:)=round(X(:,index),6); % 这是所有点的x坐标
Dataxy(2:2:end,:)=round(Y(:,index),6); %             y坐标
Datav=round(V(:,index),6); % 这些记录所有点的速度信息

% 写入文件
filename = 'result1_test.xlsx';

% 创建表头
headers = cell(1, length(index));
for i = 1:length(index)
    headers{i} = sprintf('%.1fs', index(i) - 1);
end

% 创建行名
N = 223; % 总节数
row_names = cell(2*N + 2, 1);
row_names{1} = '龙头x (m)';
row_names{2} = '龙头y (m)';
for i = 2:N-1
    row_names{2*i-1} = sprintf('第%d节龙身x (m)', i-1);
    row_names{2*i} = sprintf('第%d节龙身y (m)', i-1);
end
row_names{2*N-1} = '龙尾x (m)';
row_names{2*N} = '龙尾y (m)';
row_names{2*N+1} = '龙尾（后）x (m)';
row_names{2*N+2} = '龙尾（后）y (m)';

% 写入表头
writecell(headers, filename, 'Sheet', 1, 'Range', 'B1');
writecell(row_names, filename, 'Sheet', 1, 'Range', 'A2');

% 写入位置数据
Dataxy = round(Dataxy, 6);
writematrix(Dataxy, filename, 'Sheet', 1, 'Range', 'B2');

% 写入速度数据
writecell(headers, filename, 'Sheet', 2, 'Range', 'B1');
writecell(row_names, filename, 'Sheet', 2, 'Range', 'A2');
Datav = round(Datav, 6);
writematrix(Datav, filename, 'Sheet', 2, 'Range', 'B2');




%% 下面再看论文中的格式怎么导出
nn2=60/dt;
index2=1:nn2:length(tt); % 找出特定的时间点,来记录
index_row=[1 2:50:224 224];
Dataxy2=zeros(2*length(index_row),length(index2)); % 这些记录所有点的位置信息
Dataxy2(1:2:end,:)=round(X(index_row,index2),6) % 这是所有点的x坐标
Dataxy2(2:2:end,:)=round(Y(index_row,index2),6) %             y坐标
Datav2=round(V(index_row,index2),6) % 这些记录所有点的速度信息

%%%%%%%%% 以上全部计算完成！%%%%%%%%%%%%

%% 定义子函数：用来求解螺线上任意一个孔所在点坐标知道的情况下，如何求下一个孔所在位置
function theta=solve_theta(luoju,x1,y1,theta1,d) % 知道等距螺线上任意一点x1,y1,求在这同一条螺线上且与(x1,y1)这点相距为d的点的角度,要求这个角度一定要大于(x1,y1)这点的角度theta1
k=luoju/2/pi;
fun=@(theta)(k*theta.*cos(theta)-x1).^2+(k*theta.*sin(theta)-y1).^2-d^2;% 利用距离
q=0.01;
options = optimoptions('fsolve','Display','off'); % 不提示结果
theta=fsolve(fun,theta1+q,options); % 以前一个孔对应的角度值为基准进行非线性方程求零点，为什么加0.1？很简单，不想找比它还小的角度
while theta<=theta1 || abs(k*theta-k*theta1)>luoju/2 % 如果求解得到的角度比前面的孔对应的角度还小，或者新求出的孔与前一个孔不在一条螺线（在其他螺线）
    q=q+0.1;
    theta=fsolve(fun,theta+q,options); % 重新求一个角度 
end % 直到满足：这个点与前面孔在一条螺线，且这个点角度比前面孔的角度大!
end

