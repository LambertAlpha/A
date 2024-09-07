%% VX公众号：Matlab techniques出品，谨防假冒！
clc;close all
warning off
luoju=1.7; % 螺距
k=luoju/2/pi; % 螺线方程的系数 r=k theta
L1=341e-2;
D1=L1-27.5e-2*2; % 龙头把手两个孔之间的距离
L2=220e-2;
D2=L2-27.5e-2*2; % 其他凳子把手两个孔之间的距离
v0=1; % 头节点速度

%% 盘出螺线与盘入螺线中心对称：
% 先画出部分螺线
theta=5*2*pi:-0.01:0*pi;
r=k*theta;
x=r.*cos(theta);
y=r.*sin(theta);
figure(1)
set(gcf,'Position',[200 200 600 600]);
c1=rand(1,3);
plot(x,y,'-','Color',c1,'LineWidth',1.3)
axis equal
grid on
xlabel('x')
ylabel('y')
set(gca,'FontSize',18) % 改变坐标轴字体大小好看一点
hold on

theta=theta-pi;
r2=k*(theta+pi); % 旋转180°
x2=r2.*cos(theta);
y2=r2.*sin(theta);
c2=rand(1,3);
plot(x2,y2,'m','Color',c2,'LineWidth',1.3)
% 绘制掉头区域
R=4.5; % 掉头区域半径
x_diao=R*cos(theta);
y_diao=R*sin(theta);
c3=sort(rand(1,3));
plot(x_diao,y_diao,'Color',c3,'LineWidth',2)
% legend('盘入螺线','盘出螺线','调头边界')

%% 写出入螺线和出螺线与圆形边界的交点角度/交点处的切线斜率等几何信息
theta_ru=R/k;
theta_chu=R/k-pi;   % 出螺线是：r=k(theta+pi)=R,所以theta_chu=R/k-pi
slope=(k*sin(theta_ru)+R*cos(theta_ru))/(k*cos(theta_ru)-R*sin(theta_ru)); % 入射点和出射点处的斜率是一样的（因为中心对称！）
theta_max_1=atan(-1/slope)+pi; % 切点到圆C1中心连线的角度（大于pi VX公众号：Matlab techniques出品）
theta_dengyao=atan(tan(theta_ru))+pi-theta_max_1; % 这是两个等腰三角形的底角
rC1_C2=R/cos(theta_dengyao); % 几何关系，得到r1+r2的值
rC2=rC1_C2/3;rC1=rC2*2; % C1半径是C2的2倍;
phi=2*theta_dengyao; % C1圆周到半圆的亏损(当然也是C2的亏损)
SC1=rC1*(pi-phi); SC2=rC2*(pi-phi); % 两段圆弧的弧长，都是定值！
theta_min_1=theta_max_1-SC1/rC1;  % 圆弧C1最右边界对应的角度
theta_min_2=theta_min_1-pi; % 圆弧C2最左边边界对应的角度(指向左下角)
theta_max_2=theta_min_2+SC2/rC2; % 圆弧C2最上端边界半径所指的角度值
x_C1=R*cos(theta_ru)+rC1*cos(theta_max_1-pi);
y_C1=R*sin(theta_ru)+rC1*sin(theta_max_1-pi); % 得到圆C1的圆心坐标

x_C2=R*cos(theta_chu)-rC2*cos(theta_max_2);
y_C2=R*sin(theta_chu)-rC2*sin(theta_max_2); % 得到圆C2的圆心坐标
figure(1)
hold on
plot(x_C1+rC1*cos(linspace(theta_min_1,theta_max_1,50)),y_C1+rC1*sin(linspace(theta_min_1,theta_max_1,50)),'r','LineWidth',2)
plot(x_C1,y_C1,'r*')
plot(x_C2+rC2*cos(linspace(theta_min_2,theta_max_2,50)),y_C2+rC2*sin(linspace(theta_min_2,theta_max_2,50)),'b','LineWidth',2)
plot(x_C2,y_C2,'b*')
axis equal
% 以上给出了两段切圆弧的所有几何信息
%% 盘入曲线上头节点的位置求解
mydtheta=@(t,theta)1./(k*sqrt(1+theta.^2));
theta0=theta_ru; % 初始在入射点,下面考虑头部节点“反方向”盘入，也就是求-100到0这段时间头部节点的位置信息，
dt=0.1; % 时间步长
tspan=0:dt:100; % 求解时间点
[tt,theta]=ode45(mydtheta,tspan,theta0); % 龙格库塔法求解
X1=k*theta.*cos(theta);
Y1=k*theta.*sin(theta);
% figure(2)
% for i=1:10:length(theta)
%     title({['t=',num2str(tt(i))],'VX公众号Matlab techniques出品','头部第一个把手中心的轨迹(0到-100s)'})
%     plot(X1(i),Y1(i),'r.','MarkerSize',10)
%     hold on
%     axis equal
%     drawnow
% end
tt_ru=tt(end:-1:1);
XX=zeros(224,200/dt+1); % 每一行记录一个把手的位置关于时间的变化数据
YY=zeros(224,200/dt+1); % 每一行记录一个把手的位置关于时间的变化数据
TH=zeros(224,200/dt+1);
XX(1,1:length(X1))=X1(end:-1:1); % 先把入射线上头部节点关于时间的位置信息填好
YY(1,1:length(Y1))=Y1(end:-1:1); % 先把入射线上头部节点关于时间的位置信息填好
TH(1,1:length(theta))=theta(end:-1:1); % 想一下为什么要“翻转”矩阵，因为我们需要的是从-100到0，而求得的结果是0到100
% 圆弧C1上头节点的位置信息
tt_c1=dt:dt:SC1; % 在圆弧C1上,由于速度为1，所以弧长等于时间，所以弧长是这段时间的上边界
theta_C1=-tt_c1/rC1+theta_max_1;
TH(1,length(theta)+(1:length(tt_c1)))=theta_C1;
XX(1,length(X1)+(1:length(tt_c1)))=rC1*cos(theta_C1)+x_C1; % 把圆弧C1上头节点关于时间的位置信息填好
YY(1,length(Y1)+(1:length(tt_c1)))=rC1*sin(theta_C1)+y_C1; % 把圆弧C1上头节点关于时间的位置信息填好
% 圆弧C2
tt_c2=tt_c1(end)+dt:dt:SC1+SC2; % 在圆弧C2上,由于速度为1，所以弧长等于时间，所以弧长是这段时间的上边界
theta_C2=(tt_c2-SC1)/rC2+theta_min_1-pi; %
TH(1,length(theta)+length(theta_C1)+(1:length(tt_c2)))=theta_C2;
XX(1,length(X1)+length(theta_C1)+(1:length(tt_c2)))=rC2*cos(theta_C2)+x_C2; % 把圆弧C2上头节点关于时间的位置信息填好
YY(1,length(Y1)+length(theta_C1)+(1:length(tt_c2)))=rC2*sin(theta_C2)+y_C2; % 把圆弧C2上头节点关于时间的位置信息填好
% 盘出螺线上头节点的位置求解
mydtheta=@(t,theta)1./(k*sqrt(1+(theta+pi).^2)); % 注意出射的螺线方程，+pi
theta0=theta_chu; % 初始在出射点,下面考虑头部节点盘出，也就是t快到100这段时间头部节点的位置信息，
tspan=tt_c2(end)+dt:dt:100; % 求解时间点
[tt,theta2]=ode45(mydtheta,tspan,theta0); % 龙格库塔法求解
X2=k*(theta2+pi).*cos(theta2);
Y2=k*(theta2+pi).*sin(theta2);
TH(1,length(theta)+length(theta_C1)+length(tt_c2)+1:end)=theta2;
XX(1,length(theta)+length(theta_C1)+length(tt_c2)+1:end)=X2;
YY(1,length(theta)+length(theta_C1)+length(tt_c2)+1:end)=Y2; % 把盘出螺线上头节点的位置信息存储
figure(3)
set(gcf,'Position',[300 300 600 600])
clf
% hold on
% plot(XX(1,:),YY(1,:),'k-','LineWidth',1.2,'Marker','o',VX公众号：Matlab techniques出品防止假冒！'MarkerSize',6,'MarkerFaceColor','r')
% axis equal
for i=1:3:length(TH(1,:))
    title({['t=',num2str((i-1)*dt)],'VX公众号Matlab techniques出品','头部把手中心的轨迹(-100到100s)'})
    plot(XX(1,i),YY(1,i),'Marker','o','MarkerSize',3,'MarkerFaceColor','r')
    hold on
    axis equal
    axis([-10 10 -10 10])
    grid on
    drawnow
end
% 上面完成了：-100到100s内，每一个时间点下头把手的位置信息
%% 下面要完成的事情：每一个时间点处，既然头节点位置信息知道了，那么其余节点是不是可以循环求解了？
hwait=waitbar(0,'计算开始...');
t_total=-100:dt:100;
for j=1:length(t_total)
    if t_total(j)<=0 % 头节点将始终位于盘入螺线
        for i=2:N+1 % 在每一个时间点下，对每一行循环计算，意思是求出此时各个把手孔的位置信息
            d=D1*(i<=2)+D2*(i>2); % 分辨下是第一个凳子还是其他凳子，孔之间的距离不一样！
            thetaij=solve_theta1(luoju,XX(i-1,j),YY(i-1,j),TH(i-1,j),d); % 子函数求解下一个孔的角度值
            TH(i,j)=thetaij;
            XX(i,j)=k*thetaij*cos(thetaij);
            YY(i,j)=k*thetaij*sin(thetaij);
            %         waitbar(((j-1)*N+i)/(length(t_total)*N),hwait,'已经完成...')
        end
    elseif t_total(j)>0 && t_total(j)<=SC1 % 头节点始终位于圆弧C1段
        flag=2;
        for i=2:N+1
            d=D1*(i<=2)+D2*(i>2); % 分辨下是第一个凳子还是其他凳子，孔之间的距离不一样！
            if flag==2 % 仍然在圆弧1区域
                [xi,yi,thetai,flag]=solve_point_2_1(luoju,XX(i-1,j),YY(i-1,j),TH(i-1,j),d,rC1,x_C1,y_C1,theta_max_1);
                TH(i,j)=thetai;
                XX(i,j)=xi;
                YY(i,j)=yi;
            else  % 说明此时求得的点xi,yi已经过度到盘入螺线上
                thetaij=solve_theta1(luoju,XX(i-1,j),YY(i-1,j),TH(i-1,j),d); % 子函数求解下一个孔的角度值
                TH(i,j)=thetaij;
                XX(i,j)=k*thetaij*cos(thetaij);
                YY(i,j)=k*thetaij*sin(thetaij);
            end
        end
    elseif t_total(j)>SC1 && t_total(j)<=SC1+SC2 % 头节点始终位于圆弧C2段
        flag=3;
        for i=2:N+1
            d=D1*(i<=2)+D2*(i>2); % 分辨下是第一个凳子还是其他凳子，孔之间的距离不一样！
            if flag==3 % 仍然在圆弧2区域
                [xi,yi,thetai,flag]=solve_point_3_2(XX(i-1,j),YY(i-1,j),TH(i-1,j),d,rC1,x_C1,y_C1,rC2,x_C2,y_C2,theta_min_2);
                TH(i,j)=thetai;
                XX(i,j)=xi;
                YY(i,j)=yi;
            elseif flag==2 % 说明此时求得的点已经过度到圆弧C1上了
                [xi,yi,thetai,flag]=solve_point_2_1(luoju,XX(i-1,j),YY(i-1,j),TH(i-1,j),d,rC1,x_C1,y_C1,theta_max_1);
                TH(i,j)=thetai;
                XX(i,j)=xi;
                YY(i,j)=yi;
            else  % 说明此时求得的点xi,yi已经过度到盘入螺线上
                thetaij=solve_theta1(luoju,XX(i-1,j),YY(i-1,j),TH(i-1,j),d); % 子函数求解下一个孔的角度值
                TH(i,j)=thetaij;
                XX(i,j)=k*thetaij*cos(thetaij);
                YY(i,j)=k*thetaij*sin(thetaij);
            end
        end
    else  % 最后一段，头节点位于盘出螺线
        flag=4;
        for i=2:N+1
            d=D1*(i<=2)+D2*(i>2); % 分辨下是第一个凳子还是其他凳子，孔之间的距离不一样！
            if flag==4 % 在盘出螺线区域
                [xi,yi,thetai,flag]=solve_point_4_3(luoju,XX(i-1,j),YY(i-1,j),TH(i-1,j),d,rC2,x_C2,y_C2,theta_max_2);
                TH(i,j)=thetai;
                XX(i,j)=xi;
                YY(i,j)=yi;
            elseif flag==3 % 跑圆弧2区域了
                [xi,yi,thetai,flag]=solve_point_3_2(XX(i-1,j),YY(i-1,j),TH(i-1,j),d,rC1,x_C1,y_C1,rC2,x_C2,y_C2,theta_min_2);
                TH(i,j)=thetai;
                XX(i,j)=xi;
                YY(i,j)=yi;
            elseif flag==2 % 说明此时求得的点已经过度到圆弧C1上了
                [xi,yi,thetai,flag]=solve_point_2_1(luoju,XX(i-1,j),YY(i-1,j),TH(i-1,j),d,rC1,x_C1,y_C1,theta_max_1);
                TH(i,j)=thetai;
                XX(i,j)=xi;
                YY(i,j)=yi;
            else  % 说明此时求得的点xi,yi已经过度到盘入螺线上
                thetaij=solve_theta1(luoju,XX(i-1,j),YY(i-1,j),TH(i-1,j),d); % 子函数求解下一个孔的角度值
                TH(i,j)=thetaij;
                XX(i,j)=k*thetaij*cos(thetaij);
                YY(i,j)=k*thetaij*sin(thetaij);
            end
        end
    end
    waitbar(j/length(t_total),hwait,'已经完成...')
end

%% 可视化调头过程！
figure(100)
clf;
set(gcf,'Position',[200 200 600 600])
for j=1:size(XX,2)
    plot(XX(:,j),YY(:,j),'k-','LineWidth',1.2,'Marker','o','MarkerSize',6,'MarkerFaceColor','r')
    title({['t=',num2str(dt*(j-1)-100)],'VX公众号Matlab techniques出品','盘入-掉头-盘出的轨迹'})
    axis equal
    grid on
    xlabel('x')
    ylabel('y')
    axis([-15 15 -15 15])
    drawnow
%     hold off
end
% %%
% figure
% for i=1:size(XX,2)
%     plot(XX(1,i),YY(1,i),'Marker','o','MarkerSize',3,'MarkerFaceColor','r')
%     hold on
%     axis equal
%     axis([-10 10 -10 10])
%     grid on
%     drawnow
% end
%%%%%%公众号Matlab techniques出品！其他出处皆为抄袭！%%%%%%%%%%%%%%

%% 定义子函数：用来求解盘入螺线上任意一个孔所在点坐标知道的情况下，如何求下一个孔所在位置
function theta=solve_theta1(luoju,x1,y1,theta1,d) % 知道等距盘入螺线上任意一点x1,y1,求在这同一条螺线上且与(x1,y1)这点相距为d的点的角度,要求这个角度一定要大于(x1,y1)这点的角度theta1
k=luoju/2/pi;
fun=@(theta)(k*theta.*cos(theta)-x1).^2+(k*theta.*sin(theta)-y1).^2-d^2;% 利用距离约束
q=0.01;
options = optimoptions('fsolve','Display','off'); % 不提示结果
theta=fsolve(fun,theta1+q,options); % 以前一个孔对应的角度值为基准进行非线性方程求零点，为什么加0.01？很简单，不想找比它还小的角度
while theta<=theta1 || abs(k*theta-k*theta1)>luoju/2 % 如果求解得到的角度比前面的孔对应的角度还小，或者新求出的孔与前一个孔不在一条螺线（在其他螺线）
    q=q+0.1;
    theta=fsolve(fun,theta+q,options); % 重新求一个角度
end % 直到满足：这个点与前面孔在一条螺线，且这个点角度比前面孔的角度大! %%%%%%公众号Matlab techniques出品！其他出处皆为抄袭！%%%%%%%%%%%%%%
end
%% 定义子函数：用来求解盘出螺线上任意一个孔所在点坐标知道的情况下，如何求下一个孔所在位置
function theta=solve_theta2(luoju,x1,y1,theta1,d) % 知道等距盘入螺线上任意一点x1,y1,求在这同一条螺线上且与(x1,y1)这点相距为d的点的角度,要求这个角度一定要大于(x1,y1)这点的角度theta1
k=luoju/2/pi;
fun=@(theta)(k*(theta+pi).*cos(theta)-x1).^2+(k*(theta+pi).*sin(theta)-y1).^2-d^2;% 利用距离约束,注意此时关系式是盘出螺线r=k(theta+pi)
q=-0.1;
options = optimoptions('fsolve','Display','off'); % 不提示结果 %%%%%%公众号Matlab techniques出品！其他出处皆为抄袭！%%%%%%%%%%%%%%
theta=fsolve(fun,theta1+q,options); % 以前一个孔对应的角度值为基准进行非线性方程求零点，为什么加-0.01？很简单，不想找比它还大的角度(盘出螺线的前面的点角度大,后面的点角度小
while theta>=theta1 || abs(k*theta-k*theta1)>luoju/2 % 如果求解得到的角度比前面的孔对应的角度还大，或者新求出的孔与前一个孔不在一条螺线（在其他螺线）
    q=q-0.1;
    theta=fsolve(fun,theta+q,options); % 重新求一个角度
end % 直到满足：这个点与前面孔在一条螺线，且这个点角度比前面孔的角度大!
end


function [x,y,theta,flag]=solve_point_2_1(luoju,x1,y1,theta1,d,rC1,x_c,y_c,theta_max)
%% 这个子函数是：知道一个点在圆弧C1上，求它后面那个把手点在哪里，要考虑圆弧C1与盘入螺线过渡的区域！x_c y_c是圆弧的圆心坐标
k=luoju/2/pi;
delta_theta=2*asin(d/2/rC1); %%%%%%公众号Matlab techniques出品！其他出处皆为抄袭！%%%%%%%%%%%%%%
if delta_theta<=theta_max-theta1 % 这说明“铺”一块长度d的板后，后面把手点依然在圆弧C1上
    flag=2; % 此时，我们规定flag=2，表示下一个解还是在圆弧C1区域
    theta=theta1+delta_theta;
    x=x_c+rC1*cos(theta);
    y=y_c+rC1*sin(theta); % 写出这个把手点的坐标
else
    theta=solve_theta1(luoju,x1,y1,4.5/k,d); % 如果后面这个把手点不在圆弧C1上，那就只能在盘入螺线上！
    flag=1; % % 此时，我们规定flag=1，表示下一个解跑到盘入螺线上去了
    x=k*theta*cos(theta);
    y=k*theta*sin(theta); % 写下这个把手点坐标
end
end

function [x,y,theta,flag]=solve_point_3_2(x1,y1,theta1,d,rC1,x_c1,y_c1,rC2,x_c2,y_c2,theta_min)
%% 这个子函数是：知道一个点在圆弧C2上，求它后面那个把手点在哪里，要考虑圆弧C2与圆弧C1过渡的区域！
delta_theta=2*asin(d/2/rC2);
if delta_theta<=theta1-theta_min % 这说明“铺”一块长度d的板后，后面把手点依然在圆弧C2上
    flag=3; % 此时，我们规定flag=3，表示下一个解还是在圆弧C2区域 %公众号Matlab techniques出品！其他出处皆为抄袭
    theta=theta1-delta_theta;
    x=x_c2+rC2*cos(theta);
    y=y_c2+rC2*sin(theta); % 写出这个把手点的坐标
else
    di=sqrt((x1-x_c1)^2+(y1-y_c1)^2); % 已知把手点到圆弧C1圆心之距离
    delta_theta=acos((di^2+rC1^2-d^2)/2/di/rC1); % 余弦定理解三角形
    theta_C1_di=atan((y1-y_c1)/(x1-x_c1)); % 这是已知手把点与C1圆心的夹角
    theta=theta_C1_di+delta_theta; % 求出下一个把手点，从C2过度到了C1圆弧!
    flag=2; % % 此时，我们规定flag=2，表示下一个解跑到C1圆弧上去了
    x=x_c1+rC1*cos(theta);
    y=y_c1+rC1*sin(theta); % 写出这个把手点的坐标
end
end

function [x,y,theta,flag]=solve_point_4_3(luoju,x1,y1,theta1,d,rC2,x_c,y_c,theta_max)
%% 这个子函数是：知道一个点在圆弧盘出螺线上，求它后面那个把手点在哪里，要考虑圆弧C2与盘出螺线过渡的区域！x_c y_c是圆弧的圆心坐标
k=luoju/2/pi;
theta=solve_theta2(luoju,x1,y1,theta1,d); % 先假设后面这个把手点还在盘出螺线上，求得了角度theta
if theta>=4.5/k-pi % 这说明“铺”一块长度d的板后，后面把手点依然在盘入螺线上 %%%%%%公众号Matlab techniques出品！其他出处皆为抄袭！%%%%%%%%%%%%%%
    flag=4; % 此时，我们规定flag=4，表示下一个解还是在盘出区域
    x=k*(theta+pi)*cos(theta);
    y=k*(theta+pi)*sin(theta); % 写出这个把手点的坐标
else
    fun=@(t)(x_c+rC2*cos(theta_max-t)-x1).^2+(y_c+rC2*sin(theta_max-t)-y1).^2-d^2;% 利用距离约束,此时跨区！上一把手点在盘出螺线，未知把手点在圆弧C2
    q=0.1;
    options = optimoptions('fsolve','Display','off'); % 不提示结果
    delta_theta=fsolve(fun,theta1+q,options); % 按理说应该只有一个唯一解！
    theta=theta_max-delta_theta; % 写出此时把手点在圆弧C2上的角度
    flag=3; % % 此时，我们规定flag=3，表示下一个解跑到圆弧上去了
    x=x_c+rC2*cos(theta);
    y=y_c+rC2*sin(theta); % 写出这个把手点的坐标
end
end