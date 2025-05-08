import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.pyplot import MultipleLocator

# 读取边列表并构建图
edges = [tuple(line) for line in np.loadtxt("****")] #网络存放路径
G = nx.Graph()
G.add_edges_from(edges)
dataset = '****'  # 修改为当前数据集名称


'''
df = pd.read_excel("****")
G = nx.from_pandas_edgelist(df, 'source', 'target', create_using=nx.Graph())
dataset = '****'
'''
# 计算度序列和度分布
degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
degree_count = nx.degree_histogram(G)
print(degree_sequence)
print(degree_count)

# 准备绘图数据，只包括频率大于0的度
degrees_to_plot = []
frequencies_to_plot = []
for degree, freq in enumerate(degree_count):
    if freq > 0:
        degrees_to_plot.append(degree)
        frequencies_to_plot.append(freq)

# 绘图
fig, ax = plt.subplots(1, 1, figsize=(13, 11), dpi=260)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax.tick_params(axis='both', direction='in')
plt.title(dataset, fontproperties='Times New Roman', size=10)

# 将X轴视为分类变量
xticks_positions = range(len(degrees_to_plot))
xticklabels = [str(deg) for deg in degrees_to_plot]

# 动态标签显示逻辑
step = 6
total_ticks = len(xticks_positions)

if total_ticks <= step:
    show_indices = list(range(total_ticks))
else:
    show_indices = list(range(0, total_ticks, step))
    if show_indices[-1] != total_ticks - 1:
        show_indices.append(total_ticks - 1)

plt.xticks(
    [xticks_positions[i] for i in show_indices],
    [xticklabels[i] for i in show_indices],
    fontproperties='Times New Roman',
    size=10,
    rotation=45
)
plt.yticks(fontproperties='Times New Roman', size=10)
plt.xlabel("Degree", fontproperties='Times New Roman', size=10, labelpad=-5)
plt.ylabel("Frequency", fontproperties='Times New Roman', size=10)

# 绘制折线图
ax.plot(xticks_positions, frequencies_to_plot, marker='o', linestyle='-', color='#2c7fb8', linewidth=1, markersize=2)

# 为每个数据点添加数字标签（可选）
for x, y in zip(xticks_positions, frequencies_to_plot):
    ax.annotate(f'{y}',
                xy=(x, y),
                xytext=(0, 1.3),
                textcoords="offset points",
                ha='center', va='bottom', fontsize=6,
                fontproperties='Times New Roman',alpha=0.8)

plt.show()
