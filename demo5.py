import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams

# 设置字体，确保能够显示中文字符
rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体
rcParams['axes.unicode_minus'] = False  # 解决负号'-'显示为方块的问题

# 定义常量
initial_radius = 4.5  # 调头结束后的螺线起始半径 (m)
pitch = 1.7  # 盘出时的螺距 (m)
max_allowed_velocity = 2.0  # 各节板凳的最大允许速度 (m/s)
total_sections = 223  # 总板凳节数
head_length = 3.41  # 龙头长度 (m)
body_length = 2.20  # 龙身和龙尾长度 (m)
section_lengths = [head_length] + [body_length] * (total_sections - 1)  # 各节板凳长度

# 计算角速度
def calculate_angular_velocity(head_velocity, radius):
    return head_velocity / radius

# 计算螺线位置
def calculate_spiral_position(time, initial_radius, pitch, head_velocity):
    radius = initial_radius + pitch * time / (2 * np.pi)
    angle = time / radius
    x_pos = radius * np.cos(angle)
    y_pos = radius * np.sin(angle)
    return x_pos, y_pos, radius

# 计算板凳在螺线上各个时刻的位置
def calculate_section_position(time, section_index, path_x, path_y, angles):
    prev_x = path_x[section_index - 1]
    prev_y = path_y[section_index - 1]
    section_length = section_lengths[section_index]
    direction = angles[section_index - 1] + np.pi  # 板凳方向
    x_pos = prev_x + section_length * np.cos(direction)
    y_pos = prev_y + section_length * np.sin(direction)
    return x_pos, y_pos

# 模拟盘出螺线运动，找到最大速度
def simulate_spiral_out(head_velocity, time_steps):
    path_x = np.zeros((time_steps, total_sections))
    path_y = np.zeros((time_steps, total_sections))
    angles = np.zeros((time_steps, total_sections))  # 用于计算方向
    max_velocities = np.zeros(time_steps)

    # 模拟各个时刻的位置和速度
    for t in range(time_steps):
        path_x[t, 0], path_y[t, 0], radius = calculate_spiral_position(t, initial_radius, pitch, head_velocity)
        angles[t, 0] = np.arctan2(path_y[t, 0], path_x[t, 0])  # 计算龙头的方向
        # 计算每节板凳的位置
        for i in range(1, total_sections):
            path_x[t, i], path_y[t, i] = calculate_section_position(t, i, path_x[t], path_y[t], angles[t])
            angles[t, i] = np.arctan2(path_y[t, i], path_x[t, i])

        # 计算每节板凳的速度
        if t > 0:
            velocities = np.sqrt((path_x[t] - path_x[t-1])**2 + (path_y[t] - path_y[t-1])**2)  # 近似速度
            max_velocities[t] = np.max(velocities)

    return path_x, path_y, max_velocities

# 找到不超过2m/s的最大速度
def find_maximum_head_velocity(time_steps):
    for head_velocity in np.arange(0.5, 3.0, 0.1):  # 逐步增加龙头速度 步长为最后一个参数
        _, _, max_velocities = simulate_spiral_out(head_velocity, time_steps)
        if np.all(max_velocities <= max_allowed_velocity):
            print(f"找到最大龙头速度 head_velocity = {head_velocity:.2f} m/s")
            return head_velocity
    return None

# 模拟盘出过程,找到最大龙头速度
time_steps = 500
max_head_velocity = find_maximum_head_velocity(time_steps)

# 可视化螺线运动和速度
def plot_spiral_out(path_x, path_y, max_velocities):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))
    ax1.set_aspect('equal')
    # 绘制螺线轨迹
    for t in range(0, time_steps, 20):  # 每隔20秒绘制一次
        ax1.plot(path_x[t, :], path_y[t, :], label=f't={t}s')
    ax1.set_title('舞龙队螺线盘出路径')
    ax1.set_xlabel('x位置(m)')
    ax1.set_ylabel('y位置(m)')
    ax1.legend()

    # 绘制速度随时间的变化
    ax2.plot(np.arange(time_steps), max_velocities, label='最大速度')
    ax2.axhline(y=max_allowed_velocity, color='r', linestyle='--', label='最大允许速度2m/s')
    ax2.set_title('每节板凳的最大速度随时间变化')
    ax2.set_xlabel('时间(s)')
    ax2.set_ylabel('速度(m/s)')
    ax2.legend()
    plt.grid(True)
    plt.show()

# 绘制找到最大速度后的路径和速度
if max_head_velocity:
    path_x, path_y, max_velocities = simulate_spiral_out(max_head_velocity, time_steps)
    plot_spiral_out(path_x, path_y, max_velocities)
