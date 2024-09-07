%% VX公众号：Matlab techniques出品，谨防假冒！
clc;close all;clear
warning off
luoju_value=[50:-0.5:42]*1e-2; % 螺距取值
theta_min=zeros(size(luoju_value)); % 每一个螺距对应一个theta_min,这是头部节点所能到达的最小角度（存在一个角度，不能再往里入了因为碰撞了）
for k_luo=1:length(luoju_value)
    luoju=luoju_value(k_luo);
    k=luoju/2/pi; % 螺线方程的系数 r=k theta
    L1=341e-2;
    D1=L1-27.5e-2*2; % 龙头把手两个孔之间的距离
    L2=220e-2;
    D2=L2-27.5e-2*2; % 其他凳子把手两个孔之间的距离
    
    %% 先画出部分螺线
    % theta=16*2*pi:-0.01:0*pi;
    % r=k*theta;
    % x=r.*cos(theta);
    % y=r.*sin(theta);
    % figure(1)
    % set(gcf,'Position',[200 200 600 600]);
    % plot(x,y,'-')VX公众号：Matlab techniques出品，谨防假冒！
    % axis equal
    % grid on
    % xlabel('x')
    % ylabel('y')
    % hold on
    % plot(4.5*cos(theta),4.5*sin(theta),'y','LineWidth',2)
    %% 第一步，先确定在当前的螺距下，龙头第一个把手的位置从哪里开始,假设一个位置开始
    mydtheta=@(t,theta)-1./(k*sqrt(1+theta.^2));
    % outer_R=9/2*4/3; % 从6米的位置开始
    N0=ceil(4.5/luoju)+3; %　从4.5米圆的外面的第三条螺线开始出发
    theta0=2*pi*N0; % 初始位置时候的角度
    dt=0.2; % 时间步长
    flag=0;
    step=0;
    
    N=223; % 龙头+龙身+龙尾总的个数
    X=nan*zeros(N+1,3); % 记录每个把手点在一个时间区间内的值
    Y=nan*zeros(N+1,3); % 每一行代表每个凳子的前把手孔的位置在各个时间点处的值,因为尾部还有一个孔，所以一共223+1个，所以X和Y对应223+1行
    Theta=nan*zeros(N+1,3); % 记录每个孔在时间区间的位置对应的角度theta（这为了求速度的）
    Theta(1,3)=theta0; % 这是头把手初始时刻已知的角度值
    
    while flag==0 % 如果flag=0，说明一直没有接触，可以继续往前进
        step=step+1;
        X(:,1)=X(:,3);
        Y(:,1)=Y(:,3);
        Theta(:,1)=Theta(:,3); % 开始下一时间步之前，把当前时间区间末尾处的各个把手位移值当作下一个区间开始的数据，所以进行这样的赋值
        tspan=[0,dt/2,dt]; % 求解下一个步长下的位置点
        [tt,theta]=ode45(mydtheta,tspan,Theta(1,1)); % 龙格库塔法求解
        X1=k*theta.*cos(theta);
        Y1=k*theta.*sin(theta); % 下一个时刻，头部把手所在位置坐标
        %     for i=1:10:length(theta)
        %         title({['t=',num2str(tt(i))],'VX公众号Matlab techniques出品','头部第一个把手中心的轨迹'})
        %         plot(X1(i),Y1(i),'r.','MarkerSize',10)
        %         drawnow
        %     end
        %     hwait=waitbar(0,'计算开始...')
        %% 第二步,确定每个时间点下，头部凳子的后面一个孔（也要在螺线上）,以及龙身和龙尾凳子各个孔所在位置(都要在螺线上）
        X(1,:)=X1;
        Y(1,:)=Y1;% 第一行已知了，上面ode45求得的头部第一个把手位置数据
        Theta(1,:)=theta; % 第一行已知了，第一个把手的角度数据，上面求了
        for j=2:length(tt)
            for i=2:N+1 % 在下一个时间点下,求出此时各个把手孔的位置信息
                d=D1*(i<=2)+D2*(i>2); % 分辨下是第一个凳子还是其他凳子，孔之间的距离不一样！
                thetaij=solve_theta(luoju,X(i-1,j),Y(i-1,j),Theta(i-1,j),d); % 子函数求解下一个孔的角度值
                Theta(i,j)=thetaij;
                X(i,j)=k*thetaij*cos(thetaij);
                Y(i,j)=k*thetaij*sin(thetaij);
            end
        end
        %     hp=plot(X(:,end),Y(:,end),'k-','LineWidth',1.2,'Marker','o','MarkerSize',6,'MarkerFaceColor','r'); % 更新下此时的龙
        %     title({['t=',num2str(300+step*dt)],'VX公众号Matlab techniques出品','继续盘入，轨迹'})
        %     drawnow
        % 以上，计算完成了下一个时间点下，各个位置点 下面要判断是否接触！公众号Matlab techniques发布，其他出处皆为抄袭！
        for i=1:round(N/2)
            x_1=X(i,2);x_2=X(i+1,2);
            y_1=Y(i,2);y_2=Y(i+1,2); % 意思是，对于第一块板，先确定此时两个孔的坐标
            theta_1=Theta(i,2);
            theta_2=Theta(i+1,2);  % 分别是这两个孔的角度
            %         if isnan(theta_2) % 如果判断到了初始盘入点，剩余的点就不用判断了
            %             break;
            %         end
            index1=find((theta_1+2*pi-Theta(:,2))>0); % 从外面一层 里找危险点
            index1=index1(end-2:end); % 找到距离前把手最近的外面一层的三个点的指标
            index2=find(Theta(:,2)-(theta_2+2*pi)>0);
            if isempty(index2)
                break; % 如果到达了盘入口附近，直接不用判断了，开始下一个时刻位置更新
            else
                index2=index2(1:min(3,length(index2)));     % 找到距离后把手最近的外面一层的三个点的指标 (可能找不到三个点，因为快到初始盘入点时，上面找不到三个最近点了）
            end
            index_i=index1(1):index2(end); % 这是当前全部要考虑的外面一层的把手点的位置指标,从小往大（逆时针排布)
            n=10;m=20; % 板的均匀离散点数量,太密精度高但是运算慢
            for kk=1:length(index_i)-1 % 开始判断第i块板是否与临近的几块板接触
                X2_1=[X(index_i(kk),2);Y(index_i(kk),2)]; % 这是第index_i(kk)把手点所在坐标
                X2_2=[X(index_i(kk+1),2);Y(index_i(kk+1),2)]; % 这是第index_i(kk+1)把手点所在坐标，它两形成第二个长方形，下面就是要判断第i块板是否与这个玩意相交
                panduan=find_if_intersect(L1*(i<=1)+L2*(i>1),[x_1;y_1],[x_2;y_2],L2,X2_1,X2_2,n,m); % 调用子函数，判断是否相交，如果是，返回非空集合panduan，若没有相交，panduan为空集
                if ~isempty(panduan)
                    flag=1; % 如果任何一个有相交，设置flag为1,不用再判断了
                    break;
                end
            end
            if flag==1  % 如果判断哪一块有接触，直接不用继续下一个时刻了，这个时刻就是最大盘入时刻
                break;
            end
        end
        
        k_luo
        step
        
        %     delete(hp)
        %     close(hwait)
    end
    theta_min(k_luo)=Theta(1,end) % 将终点所在位置的角度存储
end
%% 这是把最后一个时间步下的龙的位置更新
% hp=plot(X(:,end),Y(:,end),'k-','LineWidth',1.2,'Marker','o','MarkerSize',6,'MarkerFaceColor','r'); % 更新下此时的龙
% title({['t=',num2str(300+step*dt)],'VX公众号Matlab techniques出品','继续盘入，轨迹'})
% drawnow
figure
plot(luoju_value,theta_min)
theta=9*pi./luoju_value;
hold on
plot(luoju_value,theta)
xlim([0.42,0.5])
xlabel('螺距(m)')
legend('任意一个螺距下头部把手能进入到的最小角度','圆周与螺线交界处角度和螺距的关系')
%%%%%%公众号Matlab techniques出品！其他出处皆为抄袭！%%%%%%%%%%%%%%

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

