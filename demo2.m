clc; close all; clear;
warning off;

% 定义常量
PITCH = 55e-2; % 螺距
COEFFICIENT = PITCH / (2 * pi); % 螺线方程的系数 r = k * theta
HANDLE_DIST_HEAD = 341e-2 - 2 * 27.5e-2; % 龙头把手两个孔之间的距离
HANDLE_DIST_BODY = 220e-2 - 2 * 27.5e-2; % 其他凳子把手两个孔之间的距离

% 绘制部分螺线
theta = 16 * 2 * pi : -0.01 : 0 * pi;
radius = COEFFICIENT * theta;
x = radius .* cos(theta);
y = radius .* sin(theta);
figure(1);
set(gcf, 'Position', [200, 200, 600, 600]);
plot(x, y, '--');
axis equal;
grid on;
xlabel('x');
ylabel('y');
hold on;

% 第一步：确定300秒内，龙头第一个把手的位置演化
initial_theta = 57.032076651015522; % 初始位置时的角度,我们假设从300秒开始出发，这个角度是一问里面第300秒处的头把手的角度数据！
time_step = 0.25; % 时间步长
flag = 0;
step = 0;

num_segments = 223; % 龙头+龙身+龙尾总的个数
x_coords = nan * zeros(num_segments + 1, 3); % 记录每个把手点在一个时间区间内的值
y_coords = nan * zeros(num_segments + 1, 3); % 每一行代表每个凳子的前把手孔的位置在各个时间点处的值
theta_coords = nan * zeros(num_segments + 1, 3); % 记录每个孔在时间区间的位置对应的角度theta
theta_coords(1, 3) = initial_theta; % 这是头把手初始时刻(t=300)已知的角度值

while flag == 0 % 如果flag=0，说明一直没有接触，可以继续往前进
    step = step + 1;
    x_coords(:, 1) = x_coords(:, 3);
    y_coords(:, 1) = y_coords(:, 3);
    theta_coords(:, 1) = theta_coords(:, 3); % 开始下一时间步之前，把当前时间区间末尾处的各个把手位移值当作下一个区间开始的数据，所以进行这样的赋值
    time_span = [0, time_step / 2, time_step]; % 求解下一个步长下的位置点
    [time, theta_values] = ode45(@(t, theta) -1 / (COEFFICIENT * sqrt(1 + theta^2)), time_span, theta_coords(1, 1)); % 龙格库塔法求解
    x_positions = COEFFICIENT * theta_values .* cos(theta_values);
    y_positions = COEFFICIENT * theta_values .* sin(theta_values); % 下一个时刻，头部把手所在位置坐标

    % 第二步：确定每个时间点下，头部凳子的后面一个孔（也要在螺线上）,以及龙身和龙尾凳子各个孔所在位置(都要在螺线上）
    x_coords(1, :) = x_positions;
    y_coords(1, :) = y_positions; % 第一行已知了，上面ode45求得的头部第一个把手位置数据
    theta_coords(1, :) = theta_values; % 第一行已知了，第一个把手的角度数据，上面求了
    for j = 2 : length(time)
        for i = 2 : num_segments + 1 % 在下一个时间点下,求出此时各个把手孔的位置信息
            distance = HANDLE_DIST_HEAD * (i <= 2) + HANDLE_DIST_BODY * (i > 2); % 分辨下是第一个凳子还是其他凳子，孔之间的距离不一样！
            theta_ij = solve_theta(PITCH, x_coords(i-1, j), y_coords(i-1, j), theta_coords(i-1, j), distance); % 子函数求解下一个孔的角度值
            if theta_ij > 16 * 2 * pi % 如果下一个孔的角度超出了初始位置咋办？不记录！
                theta_coords(i, j) = nan;
                x_coords(i, j) = nan;
                y_coords(i, j) = nan;
                break;
            else
                theta_coords(i, j) = theta_ij;
                x_coords(i, j) = COEFFICIENT * theta_ij * cos(theta_ij);
                y_coords(i, j) = COEFFICIENT * theta_ij * sin(theta_ij);
            end
        end
    end
    hp = plot(x_coords(:, end), y_coords(:, end), 'k-', 'LineWidth', 1.2, 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b'); % 更新下此时的龙
    title({['t = ', num2str(300 + step * time_step)], '继续盘入，轨迹'});
    drawnow;

    % 以上，计算完成了下一个时间点下，各个位置点 下面要判断是否接触
    for i = 1 : num_segments
        x_1 = x_coords(i, 2); x_2 = x_coords(i+1, 2);
        y_1 = y_coords(i, 2); y_2 = y_coords(i+1, 2); % 意思是，对于第一块板，先确定此时两个孔的坐标
        theta_1 = theta_coords(i, 2);
        theta_2 = theta_coords(i+1, 2);  % 分别是这两个孔的角度
        if isnan(theta_2) % 如果判断到了初始盘入点，剩余的点就不用判断了
            break;
        end
        index1 = find((theta_1 + 2 * pi - theta_coords(:, 2)) > 0); % 从外面一层里面找
        index1 = index1(end-2 : end); % 找到距离前把手最近的外面一层的三个点的指标
        index2 = find(theta_coords(:, 2) - (theta_2 + 2 * pi) > 0);
        if isempty(index2)
            break; % 如果到达了盘入口附近，直接不用判断了，开始下一个时刻位置更新
        else
            index2 = index2(1 : min(3, length(index2))); % 找到距离后把手最近的外面一层的三个点的指标 (可能找不到三个点，因为快到初始盘入点时，上面找不到三个最近点了）
        end
        index_i = index1(1) : index2(end); % 这是当前全部要考虑的外面一层的把手点的位置指标,从小往大（逆时针排布)
        n = 10; m = 20; % 板的均匀离散点数量,太密精度高但是运算慢
        for kk = 1 : length(index_i) - 1 % 开始判断第i块板是否与临近的几块板接触
            X2_1 = [x_coords(index_i(kk), 2); y_coords(index_i(kk), 2)]; % 这是第index_i(kk)把手点所在坐标
            X2_2 = [x_coords(index_i(kk+1), 2); y_coords(index_i(kk+1), 2)]; % 这是第index_i(kk+1)把手点所在坐标，它两形成第二个长方形，下面就是要判断第i块板是否与这个玩意相交
            panduan = find_if_intersect(HANDLE_DIST_HEAD * (i <= 1) + HANDLE_DIST_BODY * (i > 1), [x_1; y_1], [x_2; y_2], HANDLE_DIST_BODY, X2_1, X2_2, n, m); % 调用子函数，判断是否相交，如果是，返回非空集合panduan，若没有相交，panduan为空集
            if ~isempty(panduan)
                flag = 1; % 如果任何一个有相交，设置flag为1,不用再判断了
                break;
            end
        end
        if flag == 1  % 如果判断哪一块有接触，直接不用继续下一个时刻了，这个时刻就是最大盘入时刻
            break;
        end
    end

    step
    delete(hp)
end

% 这是把最后一个时间步下的龙的位置更新
hp = plot(x_coords(:, end), y_coords(:, end), 'k-', 'LineWidth', 1.2, 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b'); % 更新下此时的龙
title({['t = ', num2str(300 + step * time_step)], '继续盘入，轨迹'});
drawnow;

% 下面来求解速度数据(有各个孔的角度数据Theta，速度数据可以用数值差分方法
velocities = -COEFFICIENT * sqrt(1 + theta_coords(:, end).^2) .* (theta_coords(:, end) - theta_coords(:, end-1)) / (time_step / 2); % 用后向差分求导数,获取最后一个时间步停止时的速度数据

% 验证一下第一个孔位置的速度是不是1，也就是V矩阵的第一行的值
figure;
plot(velocities, 'b-', 'LineWidth', 1.3);
ylim([0, 1.1]);
xlabel('时间');
ylabel('头把手的速度');
title('验证数值计算得到的速度');

% 到此，我们得到了不能再盘入时刻，各个把手的坐标，存储在X(:,end)和Y(:,end)里面，速度数据在V里面，可以自己填表(参考code1那里）
% 输出到 Excel 文件
filename = 'result2_test.xlsx';

% 创建表头
headers = {'横坐标x (m)', '纵坐标y (m)', '速度 (m/s)'};
writecell(headers, filename, 'Sheet', 1, 'Range', 'B1:D1');

% 创建行名
row_names = {'龙头';};
for i = 1 : 221
    row_names{end+1} = ['第', num2str(i), '节龙身'];
end
row_names{end+1} = '龙尾';
row_names{end+1} = '龙尾（后）';

% 写入行名
writecell(row_names', filename, 'Sheet', 1, 'Range', 'A2:A225');

% 写入数据
writematrix(x_coords(:, end), filename, 'Sheet', 1, 'Range', 'B2');
writematrix(y_coords(:, end), filename, 'Sheet', 1, 'Range', 'C2');
writematrix(velocities, filename, 'Sheet', 1, 'Range', 'D2');

% 定义子函数：用来求解螺线上任意一个孔所在点坐标知道的情况下，如何求下一个孔所在位置
function theta = solve_theta(pitch, x1, y1, theta1, distance)
    coefficient = pitch / (2 * pi);
    fun = @(theta) (coefficient * theta .* cos(theta) - x1).^2 + (coefficient * theta .* sin(theta) - y1).^2 - distance^2; % 利用距离
    q = 0.01;
    options = optimoptions('fsolve', 'Display', 'off'); % 不提示结果
    theta = fsolve(fun, theta1 + q, options); % 以前一个孔对应的角度值为基准进行非线性方程求零点，为什么加0.1？很简单，不想找比它还小的角度
    while theta <= theta1 || abs(coefficient * theta - coefficient * theta1) > pitch / 2 % 如果求解得到的角度比前面的孔对应的角度还小，或者新求出的孔与前一个孔不在一条螺线（在其他螺线）
        q = q + 0.1;
        theta = fsolve(fun, theta + q, options); % 重新求一个角度
    end % 直到满足：这个点与前面孔在一条螺线，且这个点角度比前面孔的角度大!
end
