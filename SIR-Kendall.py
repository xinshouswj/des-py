import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from collections import OrderedDict
import random
from matplotlib.pyplot import MultipleLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import ConnectionPatch
# SIR参数设置
gamma = 1  # 免疫率
step = 50  # 迭代次数
m = 100  #实验次数100

markers = ['P', '>', 'x', 'X', 'o', 's', 'D','+' , 'd']
linestyles_dict =OrderedDict( [('dashed', (0, (5, 5))), ('loosely dashed',(0,(5,10))), ('dotted',(0, (1, 1))),
                    ('dashdotted', (0, (3, 5, 1, 5))),
                    ('densely dashdotdotted',(0, (3, 1, 1, 1, 1, 1))),
                    ('loosely dashdotted', (0, (3, 10, 1, 10))),('dashdotdotted', (0, (3, 5, 1, 5, 1, 5))),
                               ('loosely dotted', (0, (1, 10))), ('long dash with offset', (5, (10, 3)))])#''' ('densely dashdotdotted',(0, (3, 1, 1, 1, 1, 1))),'''
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

def yuzhi(G):
    k = 0
    k2 = 0
    for i in list(G.nodes):
        k = k + nx.degree(G, i)/ len(G.nodes)
        k2 = k2 + (nx.degree(G, i)*nx.degree(G, i)) / len(G.nodes)
    print("平均度: ", k)
    return k / (k2 - k)


def pair_is_concordant(o1, o2):
    return (o1[0] > o2[0] and o1[1] > o2[1]) or (o1[0] < o2[0] and o1[1] < o2[1])


def pair_is_disconcordant(o1, o2):
    return (o1[0] > o2[0] and o1[1] < o2[1]) or (o1[0] < o2[0] and o1[1] > o2[1])


def pair_is_xsame(o1, o2):
    return (o1[0] == o1[1])


def pair_is_ysame(o1, o2):
    return (o2[0] == o2[1])


def tau(values1, values2, get_pairs=False):
    """
    	Computes Kendall's tau rank correlation coefficient for two	lists of observations with no ties

    		Parameters:
    			values1:	list of observations
    			values2:	list of observations
    			get_pairs:	set True to return the number of concordant and discordant pairs instead of tau
        """
    concordant_pairs = 0
    discordant_pairs = 0
    zipped_vals = list(zip(values1, values2))  # 将数列转化为数组
    for i, o1 in enumerate(zipped_vals):  # 遍历数组中的元素
        for o2 in zipped_vals[i + 1:]:
            if pair_is_concordant(o1, o2):
                concordant_pairs += 1
            if pair_is_disconcordant(o1, o2):
                discordant_pairs += 1
    n = len(zipped_vals)
    print(n)
    tau = float(concordant_pairs - discordant_pairs) / (0.5 * n * (n - 1))

    if get_pairs:
        return (concordant_pairs, discordant_pairs)
    else:
        return tau

if __name__ == '__main__':
    dataset = '****' # 网络名称

    edges = [tuple(line) for line in
             np.loadtxt("****",dtype=int)] #网络存放路径
    G = nx.Graph()
    G.add_edges_from(edges)
    adj1 = nx.adjacency_matrix(G)
    adj = adj1.todense()
    nodes_n = list(G.nodes())
    nodes_n.sort()
    '''
    df = pd.read_excel("****",dtype=int) # 网络存放路径
    G = nx.from_pandas_edgelist(df, 'source', 'target', create_using=nx.Graph())
    adj1 = nx.adjacency_matrix(G)
    adj = adj1.todense()
    nodes_n = list(G.nodes())
    nodes_n.sort()
    '''
    # 数据列表
    df = pd.read_csv('****') #多方法排序列表
    metheds = df.drop(df.columns[0], axis=1)#自动略过第一列id列
    # 读取方法文件的列名
    methods_name = np.array(metheds.columns)[0: len(metheds.columns)]
    print('文件中共有以下 ' + str(len(methods_name)) + ' 种方法')
    print(methods_name)
    plt_methods = [0, 1, 2, 3, 4, 5, 6]
    yuzhi1 = yuzhi(G)
    print(yuzhi1)
    nums = 10
    f_values = np.linspace(yuzhi1, 2*yuzhi1, num=nums)#感染率变化范围
    print(f_values)
    kendall_sizes_list = []  # 放所有试验的kendall对应值
    dfSIR1 = pd.DataFrame()
    for i in range(0, nums):
        kendall_sizes = []#存放每个感染率的肯德尔
        beta = f_values[i]
        dfSIR = pd.DataFrame()
        for j in nodes_n:
            sir_list = []
            for k in range(m):
                sir_source = [j]

                sir_values = SIR_network(G, sir_source, beta, gamma, step)
                Fc = sir_values[step - 1]  # 最终的感染范围
                if Fc > 0.000001:
                    sir_list.append(Fc)


            sir = np.mean(sir_list)  # 对100实验的输出结果求均值
            dfSIR = pd.concat([dfSIR, pd.DataFrame({'SIR': [sir]})], ignore_index=True)
        dfSIR.to_csv('****', mode='a', index=False, header=False) # 存放每次模拟的SIR值
        dfSIR1 =pd.concat([dfSIR1,dfSIR],axis=1)
        sirzhi = np.array(dfSIR)#用来计算肯德尔
        for name in methods_name:#每一个方法列都和模拟出来的列计算肯德尔
            print(name)
            x = metheds[name].values.tolist()  # 读取某一列
            kendall = tau(x, sirzhi, get_pairs=False)
            kendall_sizes.append(kendall)


        kendall_sizes_list.append(kendall_sizes)
        print(kendall_sizes_list)
        df1 = pd.DataFrame(kendall_sizes_list)
        df1.to_csv('****',mode='a', index=False, header=False)#存放每次的Kendall值
    dfk = pd.DataFrame(kendall_sizes_list)
    dfk.to_csv('****')#存放所有的Kendall值
    kendall_sizes_list = np.array(kendall_sizes_list)#将list转换为矩阵，方便后面取列，按列取值画图

    dfSIR1.to_csv('****', mode='w', index=False, header=False)#将所有模拟出的值拼结后保存

    ave_tau = []#存放每个方法在不同（10个）感染率值下肯德尔值的均值，在kendall_sizes_list中是10行7列，行代表感染率值，列代表方法
    for j in range(kendall_sizes_list.shape[1]): # 列
        sum = 0
        for i in range(kendall_sizes_list.shape[0]):#行
            sum = sum + kendall_sizes_list[i][j]
        ave_tau.append(sum / nums)
    temp = pd.DataFrame(ave_tau, index=["H-index", "Cenc", "BC", "PageRank", "IMC", "KSGC","DSE"])
    temp.to_csv('****') # 存放平均Kendell

    fig, ax = plt.subplots(1, 1, figsize=(13, 11), dpi=260)
    plt.rcParams['xtick.direction'] = 'in'  # 将刻度线方向朝内
    plt.rcParams['ytick.direction'] = 'in'
    ax.tick_params(axis='both', direction='in')  # 将刻度线方向朝内
    plt.title(dataset, fontproperties='Times New Roman', size=10)  # 设置标题及字体
    plt.xticks(fontproperties='Times New Roman', size=10)  # 设置字体
    plt.yticks(fontproperties='Times New Roman', size=10)  # 设置字体
    plt.xlabel("Infection rate", fontproperties='Times New Roman', size=10, labelpad=-2.2)  # x轴标题
    plt.ylabel("Kendall's tau", fontproperties='Times New Roman', size=10)  # y轴标题:第t个时间迭代步时的感染态和免疫态的数量之和F(t)

    ax = plt.gca()  # ax为两条坐标轴的实例
    ax.xaxis.set_major_locator(MultipleLocator(0.001))  # 把x轴的刻度间隔设置为5，并存在变量里
    ax.yaxis.set_major_locator(MultipleLocator(0.05))  # 把y轴的刻度间隔设置为5，并存在变量里
    plt_labels = []  # 存储需要plt的方法名，在图例中进行设置
    # 画图
    for i in plt_methods:
        ax.plot(f_values, kendall_sizes_list[:,i], marker=markers[i],#取kendall_sizes_list的列画图
                linestyle=linestyles[i],
                linewidth=1, ms=4)
        plt_labels.append(methods_name[i])

    # 添加图例
    ax.legend(labels=plt_labels,
              prop={'family': 'Times New Roman', 'weight': 'normal', 'size': 8},
              loc=8, ncol=8, columnspacing=2,
              bbox_to_anchor=(0.49, -0.18),  # -0.2
              frameon=False)
    plt.savefig('****') # 存放最终的SIR-Kendall图片
    plt.show()