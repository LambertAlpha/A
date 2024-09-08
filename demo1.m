clc; close all;
warning off;

% 定义常量
PITCH = 55e-2; % 螺距
COEFFICIENT = PITCH / (2 * pi); % 螺线方程的系数 r = k * theta
HANDLE_DIST_HEAD = 341e-2 - 2 * 27.5e-2; % 龙头把手两个孔之间的距离
HANDLE_DIST_BODY = 220e-2 - 2 * 27.5e-2; % 其他凳子把手两个孔之间的距离

% 绘制部分螺线
theta = 20 * 2 * pi : -0.01 : 0 * pi;
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
initial_theta = 2 * pi * 16; % 初始位置时的角度
time_step = 0.05; % 时间步长
time_span = 0 : time_step : 300; % 求解时间点
[time, theta_values] = ode45(@(t, theta) -1 / (COEFFICIENT * sqrt(1 + theta^2)), time_span, initial_theta); % 龙格库塔法求解
x_positions = COEFFICIENT * theta_values .* cos(theta_values);
y_positions = COEFFICIENT * theta_values .* sin(theta_values);
for i = 1 : 10 : length(theta_values)
    title({['t = ', num2str(time(i))], '龙头前把手的运动轨迹'});
    plot(x_positions(i), y_positions(i), 'b.', 'MarkerSize', 10);
    drawnow;
end
hwait = waitbar(0, '计算开始...');

% 第二步：确定每个时间点下，头部凳子的后面一个孔（也要在螺线上）,以及龙身和龙尾凳子各个孔所在位置(都要在螺线上）
num_segments = 223; % 龙头+龙身+龙尾总的个数
x_coords = zeros(num_segments + 1, length(x_positions));
y_coords = zeros(num_segments + 1, length(x_positions)); % 每一行代表每个凳子的前把手孔的位置在各个时间点处的值
theta_coords = zeros(num_segments + 1, length(x_positions)); % 记录每个孔在各个时刻处的位置对应的角度theta
x_coords(1, :) = x_positions;
y_coords(1, :) = y_positions; % 第一行已知了，上面求得的头部第一个把手位置数据
theta_coords(1, :) = theta_values; % 第一行已知了，第一个把手的角度数据，上面求了
for j = 1 : length(time)
    for i = 2 : num_segments + 1 % 在每一个时间点下，对每一行循环计算，意思是求出此时各个把手孔的位置信息
        distance = HANDLE_DIST_HEAD * (i <= 2) + HANDLE_DIST_BODY * (i > 2); % 分辨下是第一个凳子还是其他凳子，孔之间的距离不一样！
        theta_ij = solve_theta(PITCH, x_coords(i-1, j), y_coords(i-1, j), theta_coords(i-1, j), distance); % 子函数求解下一个孔的角度值
        theta_coords(i, j) = theta_ij;
        x_coords(i, j) = COEFFICIENT * theta_ij * cos(theta_ij);
        y_coords(i, j) = COEFFICIENT * theta_ij * sin(theta_ij);
        waitbar(((j-1) * num_segments + i) / (length(time) * num_segments), hwait, sprintf('已经完成 %.2f%%...', ((j-1) * num_segments + i) / (length(time) * num_segments) * 100));
    end
end
close(hwait);

% 下面来可视化，盘入的动态图
theta = 20 * 2 * pi : -0.01 : 0 * pi;
radius = COEFFICIENT * theta;
x = radius .* cos(theta);
y = radius .* sin(theta);
set(gcf, 'Position', [200, 200, 600, 600]);
for j = 1 : length(time)
    figure(2);
    plot(x, y, '-');
    axis equal;
    grid on;
    xlabel('x');
    ylabel('y');
    hold on;
    
    plot(x_coords(:, j), y_coords(:, j), 'k-', 'LineWidth', 1.2, 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    title({['t = ', num2str(time(j))], '盘入的轨迹'});
    drawnow;
    hold off;
end

% 下面来求解速度数据(有各个孔的角度数据Theta，速度数据可以用数值差分方法
velocities = zeros(size(x_coords)); % 每一行对应一个孔的速度在不同时间节点上的值
velocities(:, 1) = -COEFFICIENT * sqrt(1 + theta_coords(:, 1).^2) .* (theta_coords(:, 2) - theta_coords(:, 1)) / time_step; % 第一个时间点，用前向差分求导数
velocities(:, end) = -COEFFICIENT * sqrt(1 + theta_coords(:, end).^2) .* (theta_coords(:, end) - theta_coords(:, end-1)) / time_step; % 最后一个时间点，用后向差分求导数
velocities(:, 2 : end-1) = -COEFFICIENT * sqrt(1 + theta_coords(:, 2 : end-1).^2) .* (theta_coords(:, 3 : end) - theta_coords(:, 1 : end-2)) / 2 / time_step; % 中间的每一个时间点，用中心差分来求导数

% 验证一下第一个孔位置的速度是不是1，也就是V矩阵的第一行的值
figure;
plot(time, velocities(1, :), 'b-', 'LineWidth', 1.3);
ylim([0, 1.1]);
xlabel('时间');
ylabel('头把手的速度');
title('验证数值计算得到的速度');

% 到此，我们已经得到了所有把手中心的位置和速度随时间的数据,下面填表
sampling_rate = 1 / time_step;
time_indices = 1 : sampling_rate : length(time); % 找出特定的时间点,来记录
position_data = zeros(2 * (num_segments + 1), length(time_indices)); % 这些记录所有点的位置信息
position_data(1 : 2 : end, :) = round(x_coords(:, time_indices), 6); % 这是所有点的x坐标
position_data(2 : 2 : end, :) = round(y_coords(:, time_indices), 6); % 这是所有点的y坐标
velocity_data = round(velocities(:, time_indices), 6); % 这些记录所有点的速度信息

% 写入文件
filename = 'result1_test.xlsx';

% 创建表头
headers = cell(1, length(time_indices));
for i = 1 : length(time_indices)
    headers{i} = sprintf('%.1fs', time_indices(i) - 1);
end

% 创建行名
row_names = cell(2 * num_segments + 2, 1);
row_names{1} = '龙头x (m)';
row_names{2} = '龙头y (m)';
for i = 2 : num_segments - 1
    row_names{2 * i - 1} = sprintf('第%d节龙身x (m)', i - 1);
    row_names{2 * i} = sprintf('第%d节龙身y (m)', i - 1);
end
row_names{2 * num_segments - 1} = '龙尾x (m)';
row_names{2 * num_segments} = '龙尾y (m)';
row_names{2 * num_segments + 1} = '龙尾（后）x (m)';
row_names{2 * num_segments + 2} = '龙尾（后）y (m)';

% 写入表头
writecell(headers, filename, 'Sheet', 1, 'Range', 'B1');
writecell(row_names, filename, 'Sheet', 1, 'Range', 'A2');

% 写入位置数据
position_data = round(position_data, 6);
writematrix(position_data, filename, 'Sheet', 1, 'Range', 'B2');

% 写入速度数据
writecell(headers, filename, 'Sheet', 2, 'Range', 'B1');
writecell(row_names, filename, 'Sheet', 2, 'Range', 'A2');
velocity_data = round(velocity_data, 6);
writematrix(velocity_data, filename, 'Sheet', 2, 'Range', 'B2');

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
