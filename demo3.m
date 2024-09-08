clc; close all; clear;
warning off;

% 定义螺距取值范围
pitch_values = [50:-0.5:42] * 1e-2; % 螺距取值
min_theta_values = zeros(size(pitch_values)); % 每一个螺距对应一个最小角度（存在一个角度，不能再往里入了因为碰撞了）

for k_pitch = 1:length(pitch_values)
    pitch = pitch_values(k_pitch);
    coefficient = pitch / (2 * pi); % 螺线方程的系数 r = k * theta
    handle_dist_head = 341e-2 - 2 * 27.5e-2; % 龙头把手两个孔之间的距离
    handle_dist_body = 220e-2 - 2 * 27.5e-2; % 其他凳子把手两个孔之间的距离
    
    % 第一步：确定在当前的螺距下，龙头第一个把手的位置从哪里开始,假设一个位置开始
    dtheta_func = @(t, theta) -1 ./ (coefficient * sqrt(1 + theta.^2));
    N0 = ceil(4.5 / pitch) + 3; % 从4.5米圆的外面的第三条螺线开始出发
    initial_theta = 2 * pi * N0; % 初始位置时的角度
    time_step = 0.2; % 时间步长
    flag = 0;
    step = 0;
    
    num_segments = 223; % 龙头+龙身+龙尾总的个数
    x_coords = nan * zeros(num_segments + 1, 3); % 记录每个把手点在一个时间区间内的值
    y_coords = nan * zeros(num_segments + 1, 3); % 每一行代表每个凳子的前把手孔的位置在各个时间点处的值
    theta_coords = nan * zeros(num_segments + 1, 3); % 记录每个孔在时间区间的位置对应的角度theta
    theta_coords(1, 3) = initial_theta; % 这是头把手初始时刻已知的角度值
    
    while flag == 0 % 如果flag=0，说明一直没有接触，可以继续往前进
        step = step + 1;
        x_coords(:, 1) = x_coords(:, 3);
        y_coords(:, 1) = y_coords(:, 3);
        theta_coords(:, 1) = theta_coords(:, 3); % 开始下一时间步之前，把当前时间区间末尾处的各个把手位移值当作下一个区间开始的数据，所以进行这样的赋值
        time_span = [0, time_step / 2, time_step]; % 求解下一个步长下的位置点
        [time, theta_values] = ode45(dtheta_func, time_span, theta_coords(1, 1)); % 龙格库塔法求解
        x_positions = coefficient * theta_values .* cos(theta_values);
        y_positions = coefficient * theta_values .* sin(theta_values); % 下一个时刻，头部把手所在位置坐标
        
        % 第二步：确定每个时间点下，头部凳子的后面一个孔（也要在螺线上）,以及龙身和龙尾凳子各个孔所在位置(都要在螺线上）
        x_coords(1, :) = x_positions;
        y_coords(1, :) = y_positions; % 第一行已知了，上面ode45求得的头部第一个把手位置数据
        theta_coords(1, :) = theta_values; % 第一行已知了，第一个把手的角度数据，上面求了
        for j = 2:length(time)
            for i = 2:num_segments + 1 % 在下一个时间点下,求出此时各个把手孔的位置信息
                distance = handle_dist_head * (i <= 2) + handle_dist_body * (i > 2); % 分辨下是第一个凳子还是其他凳子，孔之间的距离不一样！
                theta_ij = solve_theta(pitch, x_coords(i-1, j), y_coords(i-1, j), theta_coords(i-1, j), distance); % 子函数求解下一个孔的角度值
                theta_coords(i, j) = theta_ij;
                x_coords(i, j) = coefficient * theta_ij * cos(theta_ij);
                y_coords(i, j) = coefficient * theta_ij * sin(theta_ij);
            end
        end
        
        % 以上，计算完成了下一个时间点下，各个位置点 下面要判断是否接触！
        for i = 1:round(num_segments / 2)
            x_1 = x_coords(i, 2); x_2 = x_coords(i+1, 2);
            y_1 = y_coords(i, 2); y_2 = y_coords(i+1, 2); % 意思是，对于第一块板，先确定此时两个孔的坐标
            theta_1 = theta_coords(i, 2);
            theta_2 = theta_coords(i+1, 2);  % 分别是这两个孔的角度
            
            index1 = find((theta_1 + 2 * pi - theta_coords(:, 2)) > 0); % 从外面一层 里找危险点
            index1 = index1(end-2:end); % 找到距离前把手最近的外面一层的三个点的指标
            index2 = find(theta_coords(:, 2) - (theta_2 + 2 * pi) > 0);
            if isempty(index2)
                break; % 如果到达了盘入口附近，直接不用判断了，开始下一个时刻位置更新
            else
                index2 = index2(1:min(3, length(index2))); % 找到距离后把手最近的外面一层的三个点的指标 (可能找不到三个点，因为快到初始盘入点时，上面找不到三个最近点了）
            end
            index_i = index1(1):index2(end); % 这是当前全部要考虑的外面一层的把手点的位置指标,从小往大（逆时针排布)
            n = 10; m = 20; % 板的均匀离散点数量,太密精度高但是运算慢
            for kk = 1:length(index_i) - 1 % 开始判断第i块板是否与临近的几块板接触
                X2_1 = [x_coords(index_i(kk), 2); y_coords(index_i(kk), 2)]; % 这是第index_i(kk)把手点所在坐标
                X2_2 = [x_coords(index_i(kk+1), 2); y_coords(index_i(kk+1), 2)]; % 这是第index_i(kk+1)把手点所在坐标，它两形成第二个长方形，下面就是要判断第i块板是否与这个玩意相交
                panduan = find_if_intersect(handle_dist_head * (i <= 1) + handle_dist_body * (i > 1), [x_1; y_1], [x_2; y_2], handle_dist_body, X2_1, X2_2, n, m); % 调用子函数，判断是否相交，如果是，返回非空集合panduan，若没有相交，panduan为空集
                if ~isempty(panduan)
                    flag = 1; % 如果任何一个有相交，设置flag为1,不用再判断了
                    break;
                end
            end
            if flag == 1  % 如果判断哪一块有接触，直接不用继续下一个时刻了，这个时刻就是最大盘入时刻
                break;
            end
        end
        
        k_pitch
        step
    end
    min_theta_values(k_pitch) = theta_coords(1, end); % 将终点所在位置的角度存储
end

% 绘制螺距与最小角度的关系图
figure;
plot(pitch_values, min_theta_values);
theta = 9 * pi ./ pitch_values;
hold on;
plot(pitch_values, theta);
xlim([0.42, 0.5]);
xlabel('螺距(m)');
legend('任意一个螺距下头部把手能进入到的最小角度', '圆周与螺线交界处角度和螺距的关系');

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
