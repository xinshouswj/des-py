import math
import networkx as nx
import matplotlib.pyplot as plt
import csv
import pandas as pd
import numpy as np
import copy

def betweenness(G):
    dc = nx.betweenness_centrality(G)
    return dc

# 定义h-index，思路为循环遍历图中每一个节点，首先取节点的邻居节点及对应的度，再根据度进行降序排列，
# 从最大的度开始比较，i从0开始，判断i是否大于下标为i的对应度，若小，则h=i+1，若大h=i，即判断是否有i个节点的度>=i
def h_index(G):
    '''
    @description: Calculate the h-index of node n.
    @param : Graph G, node n
    @return: h-index of node n
    '''
    hindex = dict()
    for n in G.nodes():
        hi = 0
        # Define node n's neighbors.
        ns = {m: G.degree[m] for m in G[n]}
        sns = sorted(zip(ns.keys(), ns.values()), key=lambda x: x[1], reverse=True)
        for i, m in enumerate(sns):  # enumerate是一个遍历函数，i表示下标，可设置起始值，m表示sns里的值
            if i >= m[1]:  # m[1]表示下标为i的第2个元素，即度
                hi = i
                break
            hi = i + 1
        hindex[n] = hi
    return hindex

# 定义向心中心性CenC
def has_common_neighbor(G, node1, node2):
    return bool(set(G.neighbors(node1)) & set(G.neighbors(node2)))

def CenC(G):
    cenc = dict()
    c = dict()  # 节点的约束系数
    x = 0  # 表示网络规模，十进制
    ks = nx.core_number(G)  # 节点的k-shell值
    dc = nx.algorithms.centrality.degree_centrality(G)  # cenc用的是节点的度中心性而不是节点的度
    for i in range(len(G.nodes)):  # 求网络的数量级
        if (len(G.nodes)) >= (math.pow(10, i)):
            continue
        else:
            x = i
            break

    for i in list(G.nodes):  # 求节点的约束系数
        neighbor0 = list(G.neighbors(i))  # 节点i的一阶邻居节点
        ci = 0
        for j in neighbor0:
            if has_common_neighbor(G, i, j):
                sum = 1 / nx.degree(G, i)
                common = list(nx.common_neighbors(G, i, j))
                for q in list(common):
                    sum = sum + ((1 / nx.degree(G, i)) * (1 / nx.degree(G, q)))
                ci = ci + math.pow(sum, 2)
            else:
                ci = ci + math.pow((1 / nx.degree(G, i)), 2)
        c[i] = ci
    for i in list(G.nodes):  # set是节点i一二三阶邻居的集合非重复
        set = []
        for j in list(G.neighbors(i)):  # 找i的一阶邻居节点
            if j not in set:
                set.append(j)
            for k in list(G.neighbors(j)):  # 找i的二阶邻居节点
                if k not in set:
                    set.append(k)
                for l in list(G.neighbors(k)):  # 找i的三阶邻居节点
                    if l not in set:
                        set.append(l)
        set.remove(i)
        cen = 0
        for m in set:  # 计算向心中心性
            cen = cen + 4*math.pow(math.pi, 2)*(((math.exp(x*dc[i]))+(1/c[i]))/math.pow(nx.shortest_path_length(G, i, m), 2))*(ks[m] / (math.fabs(ks[i]-ks[m])+1))
        cenc[i] = cen
    return cenc

# 定义节点收缩法，收缩节点后子图的集合
# 输入为原始状态下的网络图，输出节点收缩后的图的集合
def node_shrink(Graph):  # 节点收缩法
    Graphs = []
    for i in list(Graph.nodes):  # 计算每个节点收缩后的图
        subGraph = nx.Graph()
        neighboor_nodes = [n for n in Graph.neighbors(i)]  # i节点的邻居节点
        Graph_int = copy.deepcopy(Graph)  # 深拷贝Graph
        Graph_int.remove_node(i)  # 删除节点i
        remove_edge = list([list(i) for i in Graph_int.edges])  # 删除节点i后，图中的边
        for j in range(len(remove_edge)):  # 将节点i的邻居节点都改成节点i实现节点收缩
            for k in range(len(neighboor_nodes)):
                if remove_edge[j][0] == neighboor_nodes[k]:  # 判断边是否与i相连，即判断[,]里第一个和第二个是否等于i的邻居节点
                    remove_edge[j][0] = i
                elif remove_edge[j][1] == neighboor_nodes[k]:
                    remove_edge[j][1] = i
                else:
                    continue
        subGraph.add_edges_from(remove_edge)
        subGraph.name = str(i)  # 给图赋属性，排序用
        Graphs.append(subGraph)

    # 冒泡排序,从小到大排
    for i in range(len(Graphs)):
        # Last i elements are already in place
        for j in range(0, len(Graphs) - i - 1):
            if int(float(Graphs[j].name)) > int(float(Graphs[j + 1].name)):
                Graphs[j], Graphs[j + 1] = Graphs[j + 1], Graphs[j]
    return Graphs  # 返回排序后节点收缩图的集合

# 按照收缩图，计算节点的收缩法重要性
def IMCs(graph_int, Graphs):  # 计算重要度
    njd_int = 1 / (len(list(graph_int.nodes)) * nx.average_shortest_path_length(graph_int))  # 原始网络凝聚度
    IMCs = dict()
    i = 0
    for graph in Graphs:
        i = i + 1
        njd = 1 / (len(list(graph)) * nx.average_shortest_path_length(graph))
        IMC = 1 - njd_int / njd
        IMCs[i] = IMC
    return IMCs

# 定义KSGC
def kshell(G):
    graph = G.copy()
    importance_dict = {}
    d = dict(graph.degree())  #获取节点度
    level = min(d.values())  #初始赋值为最小的度值
    while len(graph.degree):
        importance_dict[level] = []
        while True:
            level_node_list = []
            for item in graph.degree:
                if item[1] <= level:
                    level_node_list.append(item[0])
            graph.remove_nodes_from(level_node_list)
            importance_dict[level].extend(level_node_list)
            if not len(graph.degree):
                return importance_dict
            if min(graph.degree, key=lambda x: x[1])[1] > level:
                break
        level = min(graph.degree, key=lambda x: x[1])[1]
    return importance_dict

def gDegree(G):
    """
    将G.degree()的返回值变为字典
    """
    node_degrees_dict = {}
    for i in G.degree():
        node_degrees_dict[i[0]] = i[1]
    return node_degrees_dict.copy()

def get_neigbors(g, node, depth=1):
    output = {}
    layers = dict(nx.bfs_successors(g, source=node, depth_limit=depth))
    nodes = [node]
    for i in range(1, depth + 1):
        output[i] = []
        for x in nodes:
            output[i].extend(layers.get(x, []))
        nodes = output[i]
    return output

def get_ksnode(ks):
    ks_node = {}
    for k, v in ks.items():
        for i in v:
            ks_node[i] = k
    return ks_node

def KSGC(G):
    ks = kshell(G.copy())
    ks_min = min(ks)
    ks_max = max(ks)
    k = gDegree(G)
    ks_node = get_ksnode(ks)
    m_d = nx.average_shortest_path_length(G)
    score = dict()
    for i in G.nodes():
        s = 0
        neighbor = get_neigbors(G,i,int(0.5*m_d))
        for d, nodes in neighbor.items():
            for j in nodes:
                cij = math.exp((ks_node[i] - ks_node[j]) / (ks_max - ks_min))
                s += cij * ((k[i] * k[j]) / d ** 2)
        score[i] = s
    return score

# 定义PageRank算法
def pager(G):
    pag = nx.pagerank(G)
    return pag

def rneijiedian(G, i):  # 选取DSE中的局部网络
    R = 2  # R表示距离
    jiedian = list()  # 放在R=2内的节点
    for j in list(G.nodes):
        if nx.shortest_path_length(G, source=i, target=j) <= R:
            jiedian.append(j)
    jiedian.remove(i)
    print("节点", i, "的R内节点有", jiedian)
    return jiedian

def jiedianduishu(G):  # DSE中节点处于中介位置的对数
    sum = dict()
    for n in list(G.nodes):
        sum1 = 0  # 存放邻居节点连边数
        neighbor0 = list(G.neighbors(n))
        for i in range(len(neighbor0) - 1):
            for j in range(i + 1, len(neighbor0)):
                if G.has_edge(neighbor0[i], neighbor0[j]):
                    sum1 = sum1 + 1
        if nx.degree(G, n) == 1:
            sum[n] = 0
        else:
            sum[n] = ((nx.degree(G,n)*(nx.degree(G,n)-1))/2-sum1)/((nx.degree(G,n)*(nx.degree(G,n)-1))/2)
    return sum

def DSE(G):  # DSE算法
    Dji = dict()
    bian = dict()
    ks = nx.core_number(G)
    print(ks)
    for n in list(G.nodes):
        jiedian = rneijiedian(G, n)
        dji = 0
        for m in jiedian:
            if G.has_edge(n, m):
                dji = dji + nx.degree(G, m) * (1 / (1 - math.log(1 / nx.degree(G, n))))  # 以e为底
            else:
                dj = list()
                for p in nx.all_shortest_paths(G, source=n, target=m):
                    d = 0
                    for i in range(0, len(p) - 1):
                        d = d + (1 - math.log(1 / nx.degree(G, p[i])))
                    dj.append(d)
                dji = dji + nx.degree(G, m) * (1 / min(dj))
        Dji[n] = dji / len(jiedian)  # 有效距离，已经考虑了局部网络的规模
    print("D", Dji)
    for n in list(G.nodes):
        sum = 0  # 存放边的ks值
        neighbor0 = list(G.neighbors(n))
        for i in neighbor0:
            sum = sum + min(ks[i], ks[n])
        bian[n] = sum / nx.degree(G, n)
    print(bian)

    dse = dict()  # 存放每节点的方法值
    duishu = jiedianduishu(G)  # 结构洞位置的对数
    for n in G.nodes:
        zhishu = Dji[n] + duishu[n] + bian[n]
        dse[n] = zhishu
    return dse


def yuzhi(G):
    k = 0
    k2 = 0
    for i in list(G.nodes):
        k = k + nx.degree(G, i) / len(G.nodes)
        k2 = k2 + (nx.degree(G, i) * nx.degree(G, i)) / len(G.nodes)
    print("平均度: ", k)
    return k / (k2 - k)


if __name__ == "__main__":
    edges = [tuple(line) for line in np.loadtxt("****", dtype=int)] # 网络的路径
    G = nx.Graph()
    G.add_edges_from(edges)
    adj1 = nx.adjacency_matrix(G)
    adj = adj1.todense()
    pos = nx.spring_layout(G)
    # nx.draw_networkx(G,pos, node_size=500, font_size=15, node_color='red', with_labels=True)
    # plt.show()

    '''
    df = pd.read_excel("****") #网络为.xlsx格式
    G = nx.from_pandas_edgelist(df, 'source', 'target',create_using=nx.Graph())
    # nx.draw_networkx(G, node_size=400, font_size=10, node_color='white', with_labels=True)
    # plt.show()
    '''

    print("网络加载完成")
    print("节点数:", len(G.nodes))
    print("边数:", len(G.edges))
    print("平均聚集系数: ", nx.average_clustering(G))
    print("平均路径长度: ", nx.average_shortest_path_length(G))
    beta = yuzhi(G)
    print("阈值为： ", beta)

    hin = h_index(G)  # 1.局部
    print("h-index为: ", hin)
    data_with_str_keys = {str(k): v for k, v in hin.items()}  # 将每部分得到的值进行保存
    df = pd.DataFrame.from_dict(data_with_str_keys, orient='index', columns=['h-index'])
    # 重置索引
    df = df.reset_index().rename(columns={'index': 'Key'})
    df.to_excel('****', index=False) # 保存路径

    Cenc = CenC(G)  # 2.半局部
    print("向心中心性: ", Cenc)
    data_with_str_keys = {str(k): v for k, v in Cenc.items()}
    df = pd.DataFrame.from_dict(data_with_str_keys, orient='index', columns=['Cenc'])
    # 重置索引
    df = df.reset_index().rename(columns={'index': 'Key'})
    df.to_excel('****', index=False)

    be = betweenness(G)  # 3.全局
    print("介数中心性为：", be)
    data_with_str_keys = {str(k): v for k, v in be.items()}
    df = pd.DataFrame.from_dict(data_with_str_keys, orient='index', columns=['BC'])
    # 重置索引
    df = df.reset_index().rename(columns={'index': 'Key'})
    df.to_excel('****', index=False)

    par = pager(G)  # 4.随机游走
    print("pagerank值: ", par)
    data_with_str_keys = {str(k): v for k, v in par.items()}
    df = pd.DataFrame.from_dict(data_with_str_keys, orient='index', columns=['pagerank'])
    # 重置索引
    df = df.reset_index().rename(columns={'index': 'Key'})
    df.to_excel('****', index=False)

    shrink_Graphs = node_shrink(G)  # 5.连通性
    IMC = IMCs(G, shrink_Graphs)
    print("节点收缩法值: ", IMC)
    data_with_str_keys = {str(k): v for k, v in IMC.items()}
    df = pd.DataFrame.from_dict(data_with_str_keys, orient='index', columns=['IMC'])
    # 重置索引
    df = df.reset_index().rename(columns={'index': 'Key'})
    df.to_excel('****', index=False)

    ksgc = KSGC(G)  # 6.改进重力模型
    print("KSGC值: ", ksgc)
    data_with_str_keys = {str(k): v for k, v in ksgc.items()}
    df = pd.DataFrame.from_dict(data_with_str_keys, orient='index', columns=['ksgc'])
    # 重置索引
    df = df.reset_index().rename(columns={'index': 'Key'})
    df.to_excel('****', index=False)

    s = jiedianduishu(G)
    ds = DSE(G)  # 7.半局部，节点信息加连边
    print("本文方法", ds)

    data_with_str_keys = {str(k): v for k, v in ds.items()}
    df = pd.DataFrame.from_dict(data_with_str_keys, orient='index', columns=['DSE'])
    # 重置索引
    df = df.reset_index().rename(columns={'index': 'Key'})
    df.to_excel('****', index=False)

    # 将字典型结果保存为csv文件
    dict_list = [hin, Cenc, be, par, IMC, ksgc, ds]
    # 获取所有字典的键名
    keys = set().union(*(m.keys() for m in dict_list))
    # 打开 CSV 文件并将字典按行输出到文件中
    with open('****', mode='w', newline='') as file:  #所以算法计算值的保存路径
        writer = csv.writer(file)
        # 写入表头
        writer.writerow(['Id', 'H-index', 'Cenc', 'BC', 'PageRank', 'IMC', 'KSGC', 'DSE'])
        # 按行输出字典
        for key in keys:
            row = [key]
            for m in dict_list:
                if key in m:
                    row.append(m[key])
                else:
                    row.append('')
            writer.writerow(row)
