import networkx as nx
import pandas as pd
from matplotlib.pyplot import MultipleLocator
import random
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import ConnectionPatch

# SIR参数设置
beta = 0.26# 感染率
gamma = 1  # 免疫率
step = 50  # 迭代次数
n = 500  #实验次数500

markers = ['P', '>', 'x', 'X', 'o', 's', 'D','+' , 'd']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#8c564b', '#e377c2']
linestyles_dict =OrderedDict( [('dashed', (0, (5, 5))), ('loosely dashed',(0,(5,10))), ('dotted',(0, (1, 1))),
                    ('dashdotted', (0, (3, 5, 1, 5))),
('densely dashdotdotted',(0, (3, 1, 1, 1, 1, 1))),
                    ('loosely dashdotted', (0, (3, 10, 1, 10))),('dashdotdotted', (0, (3, 5, 1, 5, 1, 5))),
                               ('loosely dotted', (0, (1, 10))), ('long dash with offset', (5, (10, 3)))])#('densely dashdotdotted',(0, (3, 1, 1, 1, 1, 1))),'''
linestyles=list(linestyles_dict.values())
plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 指定默认字体
plt.rcParams['axes.unicode_minus'] = False# 解决保存图像是负号'-'显示为方块的问题

def update_node_status(graph, node, beta, gamma):
    """
    更新节点状态
    :param graph: 网络图
    :param node: 节点序数
    :param beta: 感染率
    :param gamma: 免疫率
    """
    # 如果当前节点状态为 感染者(I) 有概率gamma变为 免疫者(R)
    if graph.nodes[node]['status'] == 'I':
        p = random.random()
        if p < gamma:
            graph.nodes[node]['status'] = 'R'
    # 如果当前节点状态为 易感染者(S) 有概率beta变为 感染者(I)
    if graph.nodes[node]['status'] == 'S':
        # 获取当前节点的邻居节点
        # 无向图：G.neighbors(node)
        # 有向图：G.predecessors(node)，前驱邻居节点，即指向该节点的节点；G.successors(node)，后继邻居节点，即该节点指向的节点。
        neighbors = list(graph.neighbors(node))
        # 对当前节点的邻居节点进行遍历
        for neighbor in neighbors:
            # 邻居节点中存在 感染者(I)，则该节点有概率被感染为 感染者(I)
            if graph.nodes[neighbor]['status'] == 'I':
                p = random.random()
                if p < beta:
                    graph.nodes[node]['status'] = 'I'
                    break


def count_node(graph):
    """
    计算当前图内各个状态节点的数目
    :param graph: 输入图
    :return: 各个状态（S、I、R）的节点数目
    """
    s_num, i_num, r_num = 0, 0, 0
    for node in graph:
        if graph.nodes[node]['status'] == 'S':
            s_num += 1
        elif graph.nodes[node]['status'] == 'I':
            i_num += 1
        else:
            r_num += 1
    return s_num, i_num, r_num


def SIR_network(graph, source, beta, gamma, step):
    """
    获得感染源的节点序列的SIR感染情况
    :param graph: networkx创建的网络
    :param source: 需要被设置为感染源的节点Id所构成的序列
    :param beta: 感染率
    :param gamma: 免疫率
    :param step: 迭代次数
    """
    n = graph.number_of_nodes()  # 网络节点个数
    sir_values = []  # 存储每一次迭代后网络中感染节点数I+免疫节点数R的总和
    # 初始化节点状态
    for node in graph:
        graph.nodes[node]['status'] = 'S'  # 将所有节点的状态设置为 易感者（S）
    # 设置初始感染源
    for node in source:
        graph.nodes[node]['status'] = 'I'  # 将感染源序列中的节点设置为感染源，状态设置为 感染者（I）
    # 记录初始状态
    sir_values.append(len(source) / n)
    # 开始迭代感染
    for j in range(step):
        # 针对对每个节点进行状态更新以完成本次迭代
        for node in graph:
            update_node_status(graph, node, beta, gamma)  # 针对node号节点进行SIR过程
        s, i, r = count_node(graph)  # 得到本次迭代结束后各个状态（S、I、R）的节点数目
        sir = (i + r) / n  # 该节点的sir值为迭代结束后 感染节点数i+免疫节点数r
        sir_values.append(sir)  # 将本次迭代的sir值加入数组
    return sir_values

def ertiaotopk(G, methods, k):
    methods_name = np.array(metheds.columns)[0: len(metheds.columns)]
    df = pd.DataFrame()
    for i in range(len(methods.columns)):#i表示方法
        v0 = []
        v2 = []
        method = methods.iloc[:, i ]
        v1 = np.array(method)#对应方法节点的降序排列序列
        v0.append(v1[0])  # 将排序后的top-1节点放入top-k节点中
        v1 = np.delete(v1, 0)

        for m in range(len(v1)):  # m表示原序列下标
            for n in range(len(v0)):#n表示v0的下标
                if G.has_edge(v0[n], v1[m]):
                    v2.append(v1[m])#放邻居节点
                    break
                if (len(v0) >= k):  # k是所取top个数
                    break
                if (n>=(len(v0)-1)):
                    v0.append(v1[m])


        for j in range(len(v2)):
            if (len(v0) < k):
                v0.append(v2[j])
            else:
                break
        df1 = pd.DataFrame(v0, columns=[methods_name[i]])  # 存储各个方法的非邻居top-k节点
        df = pd.concat([df,df1],axis=1)
    return df



if __name__ == '__main__':
    dataset = '****'  # 数据集/网络名称

    edges = [tuple(line) for line in
             np.loadtxt("****",dtype=int)] # 网络存放路径
    graph = nx.Graph()
    graph.add_edges_from(edges)
    adj1 = nx.adjacency_matrix(graph)
    adj = adj1.todense()
    #nx.draw(graph, node_size=400, font_size=20, node_color='white', with_labels=True)
    #plt.show()

    '''
    df = pd.read_excel("****") # 网络存放路径
    graph = nx.from_pandas_edgelist(df, 'source', 'target', create_using=nx.Graph())
    '''

    metheds = pd.read_csv('****')  # 多种方法的排序结果
    # 读取方法文件的列名
    methods_name = np.array(metheds.columns)[0: len(metheds.columns)]
    print('文件中共有以下 ' + str(len(methods_name)) + ' 种方法')
    print(methods_name)

    # 选择需要画图显示的方法，可自行调整，也可设置为手动输入
    plt_methods = [0, 1, 2, 3, 4,5,6]

    # 可自行设置Top-K的具体k值
    #k = int(len(graph.nodes)*0.01)+1 #根据需求选取节点
    k = int(len(graph.nodes) * 0.05)+1
    #k = 1
    print(k)
    print('Top-%d...' % k)
    df = ertiaotopk(graph, metheds, k)#不取邻居节点
    sir_values_list = []  # 存储在该Top-k下各方法的sir感染情况
        # 循环所有方法
    for name in methods_name:
        sir = []
        average_sir = []
        sir_source = np.array(df[name])
        print(sir_source) #输出各方法的top-k
        for j in range(n):
            sir_source = np.array(df[name])
            sir_values = SIR_network(graph, sir_source, beta, gamma, step) # sir传播
            sir.append(sir_values) #将n次实验的结果储存为一个矩阵
        average_sir = np.mean(sir, axis=0)#求矩阵每一列的均值，求n次实验每步的均值
        sir_values_list.append(average_sir)# 存储每个方法的Sir传播情况
    temp = pd.DataFrame(sir_values_list)
    for i in plt_methods:
        print(methods_name[i], sir_values_list[i][step])#包括了初始状态

    # SIR传播曲线结果可视化
    fig, ax = plt.subplots(1, 1, figsize=(13, 11), dpi=260)
    plt.rcParams['xtick.direction'] = 'in'  # 将刻度线方向朝内
    plt.rcParams['ytick.direction'] = 'in'
    ax.tick_params( axis='both',direction='in') #将刻度线方向朝内
    plt.title(dataset + ' Top-' + str(k),fontproperties='Times New Roman',size=10)  # 设置标题及字体
    plt.xticks(fontproperties='Times New Roman', size=10)  # 设置字体
    plt.yticks(fontproperties='Times New Roman', size=10)  # 设置字体
    plt.xlabel('Time step', fontproperties='Times New Roman', size=10, labelpad=-2.5)  # x轴标题
    plt.ylabel('F(t)(%)', fontproperties='Times New Roman', size=10)  # y轴标题:第t个时间迭代步时的感染态和免疫态的数量之和F(t)

    ax = plt.gca()  # ax为两条坐标轴的实例
    ax.set_xlim(0, 51)#使得横轴刻度从零开始
    ax.xaxis.set_major_locator(MultipleLocator(5))  # 把x轴的刻度间隔设置为5，并存在变量里
    ax.yaxis.set_major_locator(MultipleLocator(5))  # 把y轴的刻度间隔设置为5，并存在变量里
    plt_labels = []  # 存储需要plt的方法名，在图例中进行设置
    # 画图
    for i in plt_methods:
        ax.plot(range(step + 1), sir_values_list[i], marker=markers[i],
                         linestyle=linestyles[i],
                         linewidth=1, ms=4)#color = colors[i],
        plt_labels.append(methods_name[i])

    # 添加图例
    ax.legend(labels=plt_labels,
                       prop={'family': 'Times New Roman','weight': 'normal', 'size': 8},
                       loc=8,ncol=8,columnspacing=2,
                       bbox_to_anchor=(0.49, -0.18),#-0.2
                       frameon = False)

    axins = inset_axes(ax, width='40%', height='38%', loc='lower right',
                               bbox_to_anchor=(0.2, 0.05, 0.8, 1),
                               bbox_transform=ax.transAxes)

    ax.tick_params(axis='both', direction='in')
    plt.xticks(fontproperties='Times New Roman', size=9)  # 设置字体
    plt.yticks(fontproperties='Times New Roman', size=9)
    axins.xaxis.set_major_locator(MultipleLocator(1))
    axins.yaxis.set_major_locator(MultipleLocator(0.01))
    for i in plt_methods:
        axins.plot(range(step + 1), sir_values_list[i], marker=markers[i],
                         linestyle=linestyles[i],
                         linewidth=1, ms=3)

    # 设置放大区间
    zone_left = 1
    zone_right = 4

    # 坐标轴的扩展比例（根据实际数据调整）
    x_ratio = 0.04  # x轴显示范围的扩展比例
    y_ratio = 0.04  # y轴显示范围的扩展比例

    # X轴的显示范围

    x = []
    for j in range(step + 1):
        x.append(j)
    xlim0 = x[zone_left] - (x[zone_right] - x[zone_left]) * x_ratio
    xlim1 = x[zone_right] + (x[zone_right] - x[zone_left]) * x_ratio

    # Y轴的显示范围
    y = []
    for i in plt_methods:
        t = np.hstack(sir_values_list[i][zone_left:zone_right+1])
        y.append(t)

    ylim0 = np.min(y) - (np.max(y) - np.min(y)) * y_ratio
    ylim1 = np.max(y) + (np.max(y) - np.min(y)) * y_ratio
    # 调整子坐标系的显示范围
    axins.set_xlim(xlim0, xlim1)
    axins.set_ylim(ylim0, ylim1)
    ax.plot([xlim0, xlim1, xlim1, xlim0, xlim0],
                    [ylim0, ylim0, ylim1, ylim1, ylim0], "black", linewidth=0.5)


    # 原图中画方框
    tx0 = xlim0
    tx1 = xlim1
    ty0 = ylim0
    ty1 = ylim1
    sx = [tx0, tx1, tx1, tx0, tx0]
    sy = [ty0, ty0, ty1, ty1, ty0]
    ax.plot(sx, sy, "black", linewidth = 0.5)

    # 画两条线
    xy = (xlim0, ylim0)
    xy2 = (xlim0, ylim1)
    con = ConnectionPatch(xyA=xy2, xyB=xy, coordsA="data", coordsB="data",
                                  axesA=axins, axesB=ax, linewidth = 0.5)
    axins.add_artist(con)

    xy = (xlim1, ylim0)
    xy2 = (xlim1, ylim1)
    con = ConnectionPatch(xyA=xy2, xyB=xy, coordsA="data", coordsB="data",
                                  axesA=axins, axesB=ax, linewidth = 0.5)
    axins.add_artist(con)

    plt.show()