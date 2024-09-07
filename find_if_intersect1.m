function flag=find_if_intersect(L1,X1_1,X1_2,L2,X2_1,X2_2,n,m)
% 意思是，知道两个长方形，分别长L1,L2,第一个长方形中两个把手点坐标为X1_1和X1_2,第二个长方形中两个把手点坐标为X2_1和X2_2
% 目的是，判断两个长方形是否有重叠！
% 公众号Matlab techniques发布，其他出处皆为抄袭！
% n m是控制离散点个数，离散的多，更精确但是耗时,n小于m,因为m是长度方向的离散，尺寸大一点
k1=(X1_1(2)-X1_2(2))/(X1_1(1)-X1_2(1)); % 第一个长方形的斜率
k1_=-1/k1; % 第一个长方形中线的垂线的斜率
k2=(X2_1(2)-X2_2(2))/(X2_1(1)-X2_2(1)); % 第二个长方形的斜率
k2_=-1/k2; % 第二个长方形中线的垂线的斜率
X1_center=(X1_1+X1_2)/2; % 第一个长方形中心点
X2_center=(X2_1+X2_2)/2; % 第二个长方形中心点
A=[k1_ -1
    k2_ -1];
P=A\[k1_*X1_center(1)-X1_center(2);k2_*X2_center(1)-X2_center(2)]; % 求解垂线的交点，新的坐标原点
vec1=X1_center-P;
vec2=X2_center-P; % 分别是两个长方形中心点与新原点连线矢量
theta1=angle(vec1(1)+1i*vec1(2)); % 
theta2=angle(vec2(1)+1i*vec2(2)); % 分别是中心点与新原点之间的夹角
delta_theta=theta1-theta2; % 两个长方形垂线之间的相对夹角！
d1=norm(vec1); % 第一个长方形中心与新原点子距离
d2=norm(vec2); % 第二个长方形中心与新原点子距离
[X,Y]=meshgrid(linspace(d1-30e-2/2,d1+30e-2/2,n),linspace(0-L1/2,0+L1/2,m)); % 均匀离散点
T=[cos(delta_theta) -sin(delta_theta)
    sin(delta_theta) cos(delta_theta)]; % 选择矩阵！
XY_new=T*[X(:) Y(:)]'; % 把所有的离散点都进行坐标旋转，转到长方形2所在的坐标系！公众号Matlab techniques发布，其他出处皆为抄袭！
flag=find(abs(XY_new(1,:)-d2)<30e-2/2 & abs(XY_new(2,:)-0)<L2/2); % 如果这些点里面有哪一个在长方形二区域里面，那么就判定相交,flag非空

%% 下面是测试这个子函数，实际运行时可以注释掉
% figure
% hold on
% plot([X1_1(1) X1_2(1)],[X1_1(2) X1_2(2)],'k')
% plot([X2_1(1) X2_2(1)],[X2_1(2) X2_2(2)],'m')
% % 公众号Matlab techniques发布，其他出处皆为抄袭！
% X1=[d1-30e-2/2 d1+30e-2/2 d1+30e-2/2 d1-30e-2/2
%     -L1/2 -L1/2 L1/2 L1/2];
% T1=[cos(theta1) -sin(theta1)
%     sin(theta1) cos(theta1)];
% X1_new=T1*X1+P;
% % figure
% fill(X1_new(1,:),X1_new(2,:),'r','FaceAlpha',0.3)
% hold on
% X2=[d2-30e-2/2 d2+30e-2/2 d2+30e-2/2 d2-30e-2/2
%     -L2/2 -L2/2 L2/2 L2/2];
% T2=[cos(theta2) -sin(theta2)
%     sin(theta2) cos(theta2)];
% X2_new=T2*X2+P;
% fill(X2_new(1,:),X2_new(2,:),'b','FaceAlpha',0.3)
% title(['flag=',num2str(flag)])
% axis equal

end
