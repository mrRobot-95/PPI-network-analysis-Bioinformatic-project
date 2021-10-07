"""
Author- Navindu Sandumina
Date_2/20/2021
Title- protein-protein interaction (PPI) network analyser
input- ppi networks files *.tsv
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

#main class input: file= ppi tsv file
class Network:
    file=''

    def __init__(self,file):

        self.file=file
        self.graph=Network.graph(self,file)

        self.protein_count=Network.pcount(self)[0]
        self.interactions = Network.pcount(self)[1]

        self.diameter=Network.diameter(self)




    #method for create protein network inputs: file= ppi.tsv file
    def graph(self,file):
        pnetwork = {}
        nodes = 0

        with open(file, 'r') as file:
            for lines in file:
                nodes += 1
                if '#' not in lines and lines != '\n':
                    edges = lines.strip().split('\t')
                    pnetwork[nodes] = edges[0:2]

        g = nx.Graph()
        for key, value in pnetwork.items():
            g.add_edge(value[0], value[1])
        return g
    #methode for count nodes and edges in network nodes are equal to nuber of proteins and edges are equal to interactions
    def pcount(self):
        g=self.graph
        protien_Count=nx.number_of_nodes(g)
        interaction=nx.number_of_edges(g)
        return protien_Count,interaction

    #methode for plot degree distribution as a histogram
    def degreePlot(self):
        g=self.graph

        degrees = sorted([d for n, d in g.degree()], reverse=True)
        plt.hist(degrees)
        plt.title("degree distribution as a histogram")
        plt.show()

        return


    #method for count diameter of given network
    def diameter(self):
        g=self.graph
        diameter=("diameter of network="+str(nx.diameter(g)))

        return diameter

    #method for plot betweenness centrality distribution histogram
    def histbetweenes(self):
        g=self.graph
        bw=nx.betweenness_centrality(g)
        plt.hist(bw)
        plt.title("betweenness centrality distribution histogram")


        return
#sub class for find details of given protein in ppi inut: file=ppi.tsv file name, protein= given protein name
class OneProtein(Network):

    def __init__(self,file,protein):
        super().__init__(file)
        self.file=file
        self.protien=protein
        self.degree=OneProtein.degree(self,protein)
        self.clusterCoefi=OneProtein.clusterCoefi(self,protein)
        self.betweenness=OneProtein.betweenCenter(self,protein)

    #methode for find the degree value of given protein
    def degree(self,protein):
        g=self.graph
        dof=g.degree(protein)
        return dof

    #methode for find the clustering coeffiecint of given protein
    def clusterCoefi(self,protein):
        g=self.graph
        undirect=g.to_undirected()
        coefi=("clustering coeffiecint of"+" "+protein+"="+str(nx.clustering(undirect,protein)))

        return coefi


    #method for find betweenness centrality of given protein
    def betweenCenter(self,protein):
        g=self.graph
        p=protein
        bw=nx.betweenness_centrality(g)
        for key,value in bw.items():
            if key == p:
                bw1=value
        bc=("betweenness centrality of"+protein+"="+str(bw1))
        return bc



#subclass for compare the two proteins of given ppi. inputs: file=ppi.tsv file, protein1 and protein2= given proteins
class TwoProteins(Network):

    def __init__(self,file,protein1,protein2):
        super().__init__(file)
        self.file=file
        self.protien1 = protein1
        self.protien2 = protein2
        self.shortest_Path=TwoProteins.shortPath(self)

    #method for find shortest path for given proteins
    def shortPath(self):
        g=self.graph
        protein1=self.protien1
        protein2=self.protien2
        path = nx.shortest_path_length(g, source=protein1 , target=protein2)
        path2 = ("Shortest path between"+" "+protein1+" " +"and"+" "+ protein2+":"+str(path))
        return path2

"""
subclass for compair two networks
input: file1 and file2 = ppi networks file name
method: compair
return: 2 barplots as average cluster coefficient and diameter of two networks and histograms of 
degree distributions and betweenness centrality distributions 
"""
class TwoNetworks(Network):
    def __init__(self,file1,file2):
        self.file1=file1
        self.file2=file2
        self.g1=Network.graph(self,file1)
        self.g2=Network.graph(self,file2)



    def compair(self):
        g1=self.g1
        g2=self.g2

        co1=nx.average_clustering(g1)
        co2=nx.average_clustering(g2)\

        nd1=nx.diameter(g1)
        nd2=nx.diameter(g2)

        x=['network1','network2']
        y1=[co1,co2]
        y2=[nd1,nd2]

        plt.figure()
        plt.subplot(1, 2, 1)
        plt.bar(x,y1)
        plt.title("average clustering coefficient of networks" )
        plt.subplot(1, 2, 2)
        plt.bar(x, y2)
        plt.title("network diameter of networks")


        bins = np.linspace(0, 100, 100)
        x3= sorted([d for n, d in g1.degree()], reverse=True)
        x4= sorted([d for n, d in g2.degree()], reverse=True)
        plt.figure()
        plt.hist(x3, bins, alpha=0.5, label='network1')
        plt.hist(x4, bins, alpha=0.5, label='network2')
        plt.legend(loc='upper right')
        plt.title("degree distributions histogram")

        bins = np.linspace(0, 100, 100)
        x3 = nx.betweenness_centrality(g1)
        x4 = nx.betweenness_centrality(g1)
        plt.figure()
        plt.hist(x3, bins, alpha=0.5, label='network1')
        plt.hist(x4, bins, alpha=0.5, label='network2')
        plt.legend(loc='upper right')
        plt.title("betweenness centrality distributions histogram of two networks")




        return




"""
main methode
obj: create object in Network class
obj2:create object in oneprotein class
obj3:create object in twoprotein class
obj4:create object in twonetwork class
uncomment for call all methodes
for the each plots first correct function was called  and then typed plt.show()
"""

if __name__ == '__main__':
    obj=Network("string_interactions.tsv")
    obj2=OneProtein("string_interactions.tsv","ERF24")
    obj3=TwoProteins("string_interactions.tsv","ARGOS","ERF24")
    obj4=TwoNetworks("string_interactions.tsv","string_interactions_2.tsv")

    print("number of proteins="+str(obj.protein_count))
    print("number of interactions="+str(obj.interactions))
    protein=obj2.protien
    print("degree of"+" "+protein+"="+str(obj2.degree))
    obj.degreePlot()
    print(obj3.shortest_Path)
    print(obj.diameter)
    print(obj2.clusterCoefi)
    print(obj2.betweenness)
    obj.histbetweenes()
    obj4.compair()
    plt.show()
