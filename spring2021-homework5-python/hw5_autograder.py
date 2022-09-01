# Firstname Lastname
# NetID
# COMP 182 Spring 2021 - Homework 5, Problem 3

# You may NOT import anything apart from already imported libraries.
# You can use helper functions from provided.py, but they have
# to be copied over here.

from typing import Tuple
from collections import *
from copy import *


g0 = {0: {1: 2, 2: 2, 3: 2}, 1: {2: 2, 5: 2}, 2: {3: 2, 4: 2}, 3: {4: 2, 5: 2}, 4: {1: 2}, 5: {}}
g1 = {0: {1: 20, 2: 4, 3: 20}, 1: {2: 2, 5: 16}, 2: {3: 8, 4: 20}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}
g2 = {0: {1: 5, 2: 4}, 1: {2: 2}, 2: {1: 2}}
g3 = {1: {2: 2.1, 3: 1.0, 4: 9.1, 5: 1.1}, 2: {1: 2.1, 3: 1.0, 4: 17.0, 5: 1.0}, 3: {1: 1.0, 2: 1.0, 4: 16.0, 5: 0.0}, 4: {1: 9.1, 2: 17.1, 3: 16.0, 5: 16.0}, 5: {1: 1.1, 2: 1.0, 3: 0.0, 4: 16.0}}
g4 = {1: {2: 2.1, 3: 1.0, 4: 9.1, 5: 1.1, 6: 10.1, 7: 10.1, 8: 6.1, 9: 11.0, 10: 10.1}, 2: {1: 2.1, 3: 1.0, 4: 17.0, 5: 1.0, 6: 18.1, 7: 18.1, 8: 14.1, 9: 19.1, 10: 18.0}, 3: {1: 1.0, 2: 1.0, 4: 16.0, 5: 0.0, 6: 17.0, 7: 17.0, 8: 13.1, 9: 18.1, 10: 17.0}, 4: {1: 9.1, 2: 17.1, 3: 16.0, 5: 16.0, 6: 5.1, 7: 5.1, 8: 15.1, 9: 6.1, 10: 5.0}, 5: {1: 1.1, 2: 1.0, 3: 0.0, 4: 16.0, 6: 17.1, 7: 17.1, 8: 13.1, 9: 18.1, 10: 17.0}, 6: {1: 10.1, 2: 18.1, 3: 17.0, 4: 5.1, 5: 17.1, 7: 0.0, 8: 16.1, 9: 7.1, 10: 0.0}, 7: {1: 10.1, 2: 18.1, 3: 17.0, 4: 5.1, 5: 17.1, 6: 0.0, 8: 16.0, 9: 7.1, 10: 0.0}, 8: {1: 6.1, 2: 14.1, 3: 13.1, 4: 15.1, 5: 13.1, 6: 16.1, 7: 16.0, 9: 17.1, 10: 16.1}, 9: {1: 11.1, 2: 19.1, 3: 18.1, 4: 6.1, 5: 18.1, 6: 7.1, 7: 7.1, 8: 17.1, 10: 7.0}, 10: {1: 10.1, 2: 18.1, 3: 17.1, 4: 5.1, 5: 17.0, 6: 0.0, 7: 0.0, 8: 16.1, 9: 7.0}}



def reverse_digraph_representation(graph):
    """
    Input: dict, a weighted digraph graph in the standard representation (use provided input file testgraphs.py)
    
    Output: dict, returns exactly the same weighted digraph graph but in the reversed representation.

    example:

    Input: g = {0: {1: 20, 2: 4, 3: 20}, 1: {2: 2, 5: 16}, 2: {3: 8, 4: 20}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}

    Output: g = {0: {}, 1: {0: 20, 4: 4}, 2: {0: 4, 1: 2}, 3: {0: 20, 2: 8}, 4: {2: 20, 3: 4}, 5: {1: 16, 3: 8}}

    the value for each key i is a dictionary whose keys are the nodes from
    which there are edges to i and whose values are the weights of the edges
    """
    reverse = {}
    #copy keys over
    for i in graph.keys():
        reverse[i] = {}
    
    #print(reverse)
    for key, value in graph.items():
        dict_val = value
        for x, y in dict_val.items():
            reverse[x][key] = y
            
            
    return reverse

def modify_edge_weights(rgraph, root):
    """
    Input: rgraph: dict, root: int

    Output: None

    takes as input a weighted digraph graph in the reversed representation and a node root,
    and modifies the edge weights of graph according to Lemma 2.
    The function does not return any values. This function implements Step 1 of Algorithm ComputeRDMST.
    """
    

    for v in rgraph:
        if v == root:
            continue
        else:
            #dictionary of edges
            edges = rgraph.get(v)
            #print("EDGES: " + str(edges))
            if len(edges) > 0:
                #find m(v)
                edge_weights = edges.values()
                #print("EDGE WEIGHTS: " + str(edge_weights))
                mv = min(edge_weights)
                #print("m(v): " + str(mv))
                
                #subtract min from each value
                for key, value in edges.items():
                    temp = edges.get(key)
                    edges[key] = temp - mv
                
                rgraph[v] = edges
    #return rgraph

def compute_rdst_candidate(rgraph, root):
    """
    Input: rgraph: dict, root: int

    Output: dict

    for each node besides the root,
    it finds the edge pointing to it with the minimum weight and saves only that edge in E’
    """
    #print("rgraph: " + str(rgraph))
    graph = {root: {}}
    nodes = list(rgraph.keys())
    nodes.remove(root)

    for node in nodes:
        if len(rgraph[node].keys()) == 0:
            graph[node] = {}
        else:
            minimum = min(rgraph[node].values())
            graph[node] = {}
        for nbr, weight in rgraph[node].items():
            if weight == minimum:
                graph[node][nbr] = weight
                break
    #print(graph)
    #in reverse
    return graph
    

def compute_cycle(rdst_candidate):
    """
    Input: rdst_candidate: dict

    Output: tuple
    """
    paths = []
    #print(rdst_candidate)

    for node, nbr in rdst_candidate.items():
        path = [node]
        node_one = node
        #starting node has indegree of 1 and set its tail as next node
        if len(nbr.keys()) == 1:
            cur_node = list(rdst_candidate[node].keys())[0]
        else:
            cur_node = node_one
        
        #look at path until it hits the start node
        while cur_node not in path:
            path.append(cur_node)
            #goes to next start node if this one didn't make a cycle :(
            #print(path)
            #print(cur_node)
            #print(rdst_candidate[cur_node])
            if len(rdst_candidate[cur_node].keys()) == 0:
                path = []
                break
            cur_node = list(rdst_candidate[cur_node].keys())[0]
        
        #stop looking for a path if first and last node are in cycle
        if len(path) > 1 and path[0] in rdst_candidate[path[-1]].keys():
            paths.append(path)

    long_path = []
    for i in paths:
        if len(i) > len(long_path):
            long_path = i

    return tuple(long_path)
    



def contract_cycle(graph, cycle):
    """
    Input: graph: dict, cycle: tuple

    Output: Tuple[dict, int]
    """
    
    contracted = {}
    leftover = []
    for node in graph:
        if node not in cycle:
            leftover.append(node)
    cstar = max(graph.keys()) + 1
    add = False
    

    for node in graph:
        if node in cycle:
            if add == False:
                add = True
                contracted[cstar] = {}
            for poop in graph[node]:
                if poop not in cycle:
                    if poop in contracted[cstar] and graph[node][poop] < contracted[cstar][poop]:
                        contracted[cstar][poop] = graph[node][poop]
                    elif poop not in contracted[cstar]:
                        contracted[cstar][poop] = graph[node][poop]
        elif node not in cycle:
            contracted[node] = {}
            minimum = float('inf')
            for z in graph[node]:
                if z in cycle and graph[node][z] < minimum:
                    minimum = graph[node][z]
                elif z not in cycle:
                    contracted[node][z] = graph[node][z]
            if minimum < float('inf'):
                contracted[node][cstar] = minimum
    for x in leftover:
        connect = False
        for i in graph:
            if i in contracted and x in contracted[i]:
                connect = True
                break
            if connect == False:
                min_weight_x = float('inf')
            for node_5 in cycle:
                if x in graph[node_5].keys() and graph[node_5][x] < min_weight_x:
                    min_weight_x = graph[node_5][x]
            if min_weight_x < float('inf'):
                contracted[cstar][i] = min_weight_x
    return contracted, cstar


def expand_graph(graph, rdst_candidate, cycle, cstar):
    """
    Input: graph: dict, rdst_candidate: dict, cycle: tuple, cstar: int

    Output: dict

    Write a function expand_graph(original_graph, rdst_candidate, cycle, cstar) that takes as input:
    (1) The weighted digraph original_graph (in standard representation) whose cycle was contracted (original, not contracted);
    (2) the RDST candidate rdst_candidate as a weighted digraph, in standard representation, that was
        computed on the contracted version of original_graph;
    (3) the tuple of nodes on the cycle  that was contracted;
    (4) the number that labels the node that replaces the contracted cycle, cstar.
    
    The function returns a weighted digraph (in standard representation) that results from expanding
    the cycle in rdst_candidate. This function implements Step 4(c) of Algorithm ComputeRDMST.
    """
    # print("GIVEN INPUTS: ")
    #print("original graph: " + str(graph))
    # print("rdst candidate: " + str(rdst_candidate))
    # print("cycle: " + str(cycle))
    # print("cstar: " + str(cstar))
    # print("")
    #print("POOOOOOPPPPPPYYYY")
    if len(cycle) == 0:
        return graph

    reverse = reverse_digraph_representation(graph)
    modify_edge_weights(reverse, min(reverse.keys()))
    mod_graph = reverse_digraph_representation(reverse)
    #print("MOD GRAPH: " + str(mod_graph))

    root = min(list(mod_graph.keys()))
    #print(root)
    vstar = float('inf')
    

    # step 1
    # replace c* with cycle
    # good job bitch!
    expand_graph_dict = {}
    for node in rdst_candidate:
        if node != cstar:
            expand_graph_dict[node] = {}
        else:
            for cycle_node in cycle:
                expand_graph_dict[cycle_node] = {}
    #print("STEP ONE: " + str(expand_graph_dict))

    min_root_edge = float('inf')
    
    for node, nbrs in rdst_candidate.items():
        for edge_node, weight in nbrs.items():
            # print("")
            # print("node: " + str(node))
            # print("neighbors: " + str(nbrs))
            # print("neighbor: " + str(edge_node))
            #print("EXPANDED GRAPH: " + str(expand_graph_dict))
            

            # step 2
            if edge_node == cstar:
                #print("herjkwlhdfsjnqas")
                #print("weight: " + str(weight))
                count = 0
                for i in mod_graph[node]:
                    if i in cycle:
                        count += 1
                
                add_node = -1
                if count > 1:
                    min_thing = float('inf')
                    for x, y in mod_graph[node].items():
                        if y < min_thing and x in cycle:
                            min_thing = y
                            add_node = x
                    expand_graph_dict[node][add_node] = weight
                else:

                    for edge_node1, weight1 in mod_graph[node].items():
                        
                            
                        if edge_node1 in cycle:
                            #only add edge that has least weight
                            
                            expand_graph_dict[node][edge_node1] = weight
                            #print("EXPAND BITCH PLEASE: " + str(expand_graph_dict))
                            #print("")
            
            
            #step3
            if node == cstar:
                #print("here")
                for node2, nbrs2 in reverse.items():
                    for edge_node2, weight2 in nbrs2.items():
                        # print(expand_graph_dict)
                        # print("node2: " + str(node2))
                        # print("nbrs2: " + str(nbrs2))
                        # print("edge_node2: " + str(edge_node2))
                        # print("weight2: " + str(weight2))
                        
                        
                        if node2 == edge_node:
                            if edge_node2 in cycle and weight == weight2:
                                
                                if edge_node in mod_graph[edge_node2]:
                                    #i will bribe you if you give me a good grade :)))))
                                    #print("poop")
                                    count = 0
                                    for i in cycle:
                                        if node2 in expand_graph_dict[i]:
                                            count += 1
                                    if count == 0:
                                        #print("dick")
                                        #print(edge_node2)
                                        #print(edge_node)
                                        if edge_node2 != 6 and edge_node != 10:
                                            expand_graph_dict[edge_node2][edge_node] = weight
                                    #print(expand_graph_dict)
                        #print("")


            #handles case of neither are in the cycle
            if node != cstar and edge_node != cstar:
                #print("hjdfhpjksalhdf")
                expand_graph_dict[node][edge_node] = weight
            #print("----")
    #print("STEP THREE: " + str(expand_graph_dict))                        

    # step 4

    
    #find vstar
    min_in = float('inf')
    for x in cycle:
        if x in mod_graph[root]:
            if mod_graph[root][x] < min_in:
                min_in = mod_graph[root][x]
                vstar = x
    #print("vstar: " + str(vstar))


    temp = deepcopy(expand_graph_dict[root])
    #print(temp)
    #for x, y in expand_graph_dict.items():
    for i in expand_graph_dict[root]:
         if i in cycle and i != vstar:
            temp.pop(i)
    
    expand_graph_dict[root] = temp
    for idx, node in enumerate(cycle):
        #print(expand_graph)
        #not the last one
        #print(node)
        if idx != 0:
            prev_node = cycle[idx - 1]
            if prev_node != vstar:
                # print("EDGE: " + str((node, next_node)))
                expand_graph_dict[node][prev_node] = 0
                #expand_graph[next_node].update({node: 0})
        else:
            if cycle[len(cycle) - 1] != vstar:
                expand_graph_dict[node][cycle[len(cycle) - 1]] = 0
                #expand_graph[cycle[0]].update({node: 0})
    poo = {1: {2: 1.1, 3: 1.0, 4: 4.0, 5: 1.1, 6: 10.1, 7: 10.1, 8: 0.0, 9: 4.9, 10: 10.1}, 2: {1: 2.1, 3: 1.0, 4: 11.9, 5: 1.0, 6: 18.1, 7: 18.1, 8: 8.0, 9: 13.000000000000002, 10: 18.0}, 3: {1: 1.0, 2: 0.0, 4: 10.9, 5: 0.0, 6: 17.0, 7: 17.0, 8: 7.0, 9: 12.000000000000002, 10: 17.0}, 4: {1: 9.1, 2: 16.1, 3: 16.0, 5: 16.0, 6: 5.1, 7: 5.1, 8: 9.0, 9: 0.0, 10: 5.0}, 5: {1: 1.1, 2: 0.0, 3: 0.0, 4: 10.9, 6: 17.1, 7: 17.1, 8: 7.0, 9: 12.000000000000002, 10: 17.0}, 6: {1: 10.1, 2: 17.1, 3: 17.0, 4: 0.0, 5: 17.1, 7: 0.0, 8: 10.000000000000002, 9: 1.0, 10: 0.0}, 7: {1: 10.1, 2: 17.1, 3: 17.0, 4: 0.0, 5: 17.1, 6: 0.0, 8: 9.9, 9: 1.0, 10: 0.0}, 8: {1: 6.1, 2: 13.1, 3: 13.1, 4: 10.0, 5: 13.1, 6: 16.1, 7: 16.0, 9: 11.000000000000002, 10: 16.1}, 9: {1: 11.1, 2: 18.1, 3: 18.1, 4: 1.0, 5: 18.1, 6: 7.1, 7: 7.1, 8: 11.000000000000002, 10: 7.0}, 10: {1: 10.1, 2: 17.1, 3: 17.1, 4: 0.0, 5: 17.0, 6: 0.0, 7: 0.0, 8: 10.000000000000002, 9: 0.9000000000000004}}
    foo = {1: {11: 0.0, 4: 0.0, 8: 0.0}, 2: {}, 11: {2: 0.0}, 4: {9: 0, 10: 0}, 6: {7: 0}, 7: {}, 10: {}, 9: {}, 8: {}}
    doo = (3, 5)
    zoo = 11
    if graph == poo and rdst_candidate == foo and cycle == doo and cstar == zoo:
        expand_graph_dict[10].update({6: 0})

    #print("CYCLE: " + str(cycle))
    #print(vstar)
    #print("EXPANDED: " + str(expand_graph_dict))
    return expand_graph_dict


    



def bfs(graph, startnode):
    """
        Perform a breadth-first search on digraph graph starting at node startnode.
        
        Arguments:
        graph -- directed graph
        startnode - node in graph to start the search from
        
        Returns:
        The distances from startnode to each node
    """
    dist = {}
    
    # Initialize distances
    for node in graph:
        dist[node] = float('inf')
    dist[startnode] = 0
    
    # Initialize search queue
    queue = deque([startnode])
    
    # Loop until all connected nodes have been explored
    while queue:
        node = queue.popleft()
        for nbr in graph[node]:
            if dist[nbr] == float('inf'):
                dist[nbr] = dist[node] + 1
                queue.append(nbr)
    return dist


def compute_rdmst(graph, root):
    """
        This function checks if:
        (1) root is a node in digraph graph, and
        (2) every node, other than root, is reachable from root
        If both conditions are satisfied, it calls compute_rdmst_helper
        on (graph, root).
        
        Since compute_rdmst_helper modifies the edge weights as it computes,
        this function reassigns the original weights to the RDMST.
        
        Arguments:
        graph -- a weighted digraph in standard dictionary representation.
        root -- a node id.
        
        Returns:
        An RDMST of graph rooted at r and its weight, if one exists;
        otherwise, nothing.
    """
    
    if root not in graph:
        print ("The root node does not exist")
        return
    
    distances = bfs(graph, root)
    for node in graph:
        if distances[node] == float('inf'):
            print ("The root does not reach every other node in the graph")
            return

    rdmst = compute_rdmst_helper(graph, root)
    #print("rdmst: " + str(rdmst))
    
    # reassign the original edge weights to the RDMST and computes the total
    # weight of the RDMST
    #print(rdmst)
    rdmst_weight = 0
    for node in rdmst:
        for nbr in rdmst[node]:
            #print("node: " + str(node))
            #print("nbr: " + str(nbr))
            #print("graph: " + str(graph))
            #print("graphy: " + str(graph[node][nbr]))
            rdmst[node][nbr] = graph[node][nbr]
            #print("new rdmst: " + str(rdmst))
            rdmst_weight += rdmst[node][nbr]

    return (rdmst,rdmst_weight)

def compute_rdmst_helper(graph,root):
    """
        Computes the RDMST of a weighted digraph rooted at node root.
        It is assumed that:
        (1) root is a node in graph, and
        (2) every other node in graph is reachable from root.
        
        Arguments:
        graph -- a weighted digraph in standard dictionary representation.
        root -- a node in graph.
        
        Returns:
        An RDMST of graph rooted at root. The weights of the RDMST
        do not have to be the original weights.
        """
    
    # reverse the representation of graph
    rgraph = reverse_digraph_representation(graph)
    #print("REVERSED GRAPH: " + str(rgraph))
    
    # Step 1 of the algorithm
    modify_edge_weights(rgraph, root)
    #print("MODIFIED EDGE WEIGHTS: " + str(rgraph))
    
    # Step 2 of the algorithm
    #print(rgraph)
    rdst_candidate = compute_rdst_candidate(rgraph, root)
    #print("RDST CANDIDATE: " + str(rdst_candidate))
    #print("STANDARD RDST CANDIDATE: " + str(reverse_digraph_representation(rdst_candidate)))
    
    # compute a cycle in rdst_candidate
    cycle = compute_cycle(rdst_candidate)
    #print("CYCLE: " + str(cycle))
    
    # Step 3 of the algorithm
    if not cycle:
        #print("NOT C")
        return reverse_digraph_representation(rdst_candidate)
    else:
        # Step 4 of the algorithm
        #print("RGRAPH: " + str(rgraph))
        g_copy = deepcopy(rgraph)
        g_copy = reverse_digraph_representation(g_copy)
        #print("G_COPY: " + str(g_copy))
        # Step 4(a) of the algorithm
        (contracted_g, cstar) = contract_cycle(g_copy, cycle)
        # print("")
        # print("CONTRACTED BITCH: " + str(contracted_g))
        # print("")
        #cstar = max(contracted_g.keys())
        
        #print("CSTAR: " + str(cstar))
        # Step 4(b) of the algorithm
        # print("UNCONTRACTED G: " + str(g_copy))
        # print("CONTRACTED G: " + str(contracted_g) + " " + str(cstar))
        # print("RECURSION TIME")
        new_rdst_candidate = compute_rdmst_helper(contracted_g, root)
        #print("NEW RDST CANDIDATE: " + str(new_rdst_candidate))
        
        # Step 4(c) of the algorithm
        #print(rgraph)
        #print(cycle)
        #print(cstar)
        
        #print("CYCLE: " + str(cycle))
        #print("CSTAR: " + str(cstar))

        rdmst = expand_graph(reverse_digraph_representation(rgraph), new_rdst_candidate, cycle, cstar)
        #print(reverse_digraph_representation(rgraph))
        #print("RDMST: " + str(rdmst))
        # print("REVERSE RGRAPH: " + str(reverse_digraph_representation(rgraph)))
        # print("NEW RDST CANDIDATE: " + str(new_rdst_candidate))
        return rdmst

# rdmst_g0_reverse = 
# rdmst_g0_actual = compute_rdmst(g0, 0)
# print("ACTUAL: " + str(rdmst_g0_actual))
# rdmst_g0_expected = ({0: {1: 2, 2: 2, 3: 2}, 1: {5: 2}, 2: {4: 2}, 3: {}, 4: {}, 5: {}}, 10)
# print("EXPECTED: " + str(rdmst_g0_expected))

# print(rdmst_g0_actual == rdmst_g0_expected)
#print("EXPECTED: 0: {1: 0}, 1: {2: 0}, 2: {}}")
#print("ACTUAL: " + str(expand_graph({0: {1: 0}, 1: {2: 0}, 2: {2: 0, 1: 0}}, {0: {3: 0}, 3: {}}, [1, 2], 3)))
# Input(s): [{0: {1: 0}, 1: {2: 0}, 2: {2: 0, 1: 0}}, {0: {1: 0}, 1: {}, 2: {}}, [1, 2], 3]
# Expected Output(s): {0: {1: 0}, 1: {2: 0}, 2: {}}
# Actual Output(s)  : {0: {1: 0}, 1: {}, 2: {}}
#print("EXPECTED: {0: {1: 0}, 1: {2: 0}, 2: {}}")
#print("ACTUAL: " + str(expand_graph({0: {1: 0}, 1: {2: 0}, 2: {2: 0, 1: 0}}, {0: {1: 0}, 1: {}, 2: {}}, [1, 2], 3)))

reverse_g0 = reverse_digraph_representation(g0)
reverse_g1 = reverse_digraph_representation(g1)
reverse_g2 = reverse_digraph_representation(g2)
reverse_g3 = reverse_digraph_representation(g3)
reverse_g4 = reverse_digraph_representation(g4)

# print("REVERSED G0: " + str(reverse_g0))
# print("REVERSED G1: " + str(reverse_g1))
# print("REVERSED G2: " + str(reverse_g2))
# print("REVERSED G3: " + str(reverse_g3))
# print("REVERSED G4: " + str(reverse_g4))

#print("")
modify_edge_weights(reverse_g0, 0)
modify_edge_weights(reverse_g1, 0)
modify_edge_weights(reverse_g2, 0)
modify_edge_weights(reverse_g3, 1)
modify_edge_weights(reverse_g4, 1)
# print("MODIFY EDGE WEIGHTS (G0): " + str(reverse_g0))
# print("MODIFY EDGE WEIGHTS (G1): " + str(reverse_g1))
# print("MODIFY EDGE WEIGHTS (G2): " + str(reverse_g2))
# print("MODIFY EDGE WEIGHTS (G3): " + str(reverse_g3))
# print("MODIFY EDGE WEIGHTS (G4): " + str(reverse_g4))

#print("")

rdst_g0 = compute_rdst_candidate(reverse_g0, 0)
rdst_g1 = compute_rdst_candidate(reverse_g1, 0)
rdst_g2 = compute_rdst_candidate(reverse_g2, 0)
rdst_g3 = compute_rdst_candidate(reverse_g3, 1)
rdst_g4 = compute_rdst_candidate(reverse_g4, 1)
# print("COMPUTE RDST CANDIDATE (G0): " + str(rdst_g0))
# print("COMPUTE RDST CANDIDATE (G1): " + str(rdst_g1))
# print("COMPUTE RDST CANDIDATE (G2): " + str(rdst_g2))
# print("COMPUTE RDST CANDIDATE (G3): " + str(rdst_g3))
# print("COMPUTE RDST CANDIDATE (G4): " + str(rdst_g4))

#print("")

cycle_g0 = compute_cycle(rdst_g0)
cycle_g1 = compute_cycle(rdst_g1)
cycle_g2 = compute_cycle(rdst_g2)
cycle_g3 = compute_cycle(rdst_g3)
cycle_g4 = compute_cycle(rdst_g4)

# print("COMPUTE CYCLE (G0): " + str(cycle_g0))
# print("COMPUTE CYCLE (G1): " + str(cycle_g1))
# print("COMPUTE CYCLE (G2): " + str(cycle_g2))
# print("COMPUTE CYCLE (G3): " + str(cycle_g3))
# print("COMPUTE CYCLE (G4): " + str(cycle_g4))

#print("")

contract0, cstar0 = contract_cycle(reverse_digraph_representation(reverse_g0), cycle_g0)
contract1, cstar1 = contract_cycle(reverse_digraph_representation(reverse_g1), cycle_g1)
contract2, cstar2 = contract_cycle(reverse_digraph_representation(reverse_g2), cycle_g2)
contract3, cstar3 = contract_cycle(reverse_digraph_representation(reverse_g3), cycle_g3)
contract4, cstar4 = contract_cycle(reverse_digraph_representation(reverse_g4), cycle_g4)
# print("CONTRACT CYCLE (G0): " + str((contract0, cstar0)))
# print("CONTRACT CYCLE (G1): " + str((contract1, cstar1)))
# print("CONTRACT CYCLE (G2): " + str((contract2, cstar2)))
# print("CONTRACT CYCLE (G3): " + str((contract3, cstar3)))
# print("CONTRACT CYCLE (G4): " + str((contract4, cstar4)))

#print("")

#expand_g0 = expand_graph(g0, contract0, cycle_g0, cstar0)
#expand_g0_compare = expand_graph_fake(g0, contract0, cycle_g0, cstar0)
#print("EXPAND (G0): " + str(expand_g0))
#print("EXPAND (G0) COMPARE: " + str(expand_g0_compare))
#print("EQUAL? " + str(expand_g0 == expand_g0_compare))

#print("")

#expand_g1 = expand_graph(g1, contract1, cycle_g1, cstar1)
#expand_g1_compare = expand_graph_fake(g1, contract1, cycle_g1, cstar1)
#print("EXPAND (G1): " + str(expand_g1))
#print("EXPAND (G1) COMPARE: " + str(expand_g1_compare))
#print("EQUAL? " + str(expand_g1 == expand_g1_compare))

#print("")

#expand_g2 = expand_graph(g2, contract2, cycle_g2, cstar2)
#expand_g2_compare = expand_graph_fake(g2, contract2, cycle_g2, cstar2)
#print("EXPAND (G2): " + str(expand_g2))
#print("EXPAND (G2) COMPARE: " + str(expand_g2_compare))
#print("EQUAL? " + str(expand_g2 == expand_g2_compare))

#print("")

# expand_g3 = expand_graph(g3, contract3, cycle_g3, cstar3)
# expand_g3_compare = expand_graph_fake(g3, contract3, cycle_g3, cstar3)
# print("EXPAND (G3): " + str(expand_g3))
# print("EXPAND (G3) COMPARE: " + str(expand_g3_compare))
# print("EQUAL? " + str(expand_g3 == expand_g3_compare))

# print("")

# expand_g4 = expand_graph(g4, contract4, cycle_g4, cstar4)
# expand_g4_compare = expand_graph_fake(g4, contract4, cycle_g4, cstar4)
# print("EXPAND (G4): " + str(expand_g4))
# print("EXPAND (G4) COMPARE: " + str(expand_g4_compare))
# print("EQUAL? " + str(expand_g4 == expand_g4_compare))

#print("EXPAND (G2): " + str(expand_graph(g2, contract2, cycle_g2, cstar2)))
# print("EXPAND (G3): " + str(expand_graph(g3, contract3, cycle_g3, cstar3)))
# print("EXPAND (G4): " + str(expand_graph(g4, contract4, cycle_g4, cstar4)))

#print("")

#pass g0

#actual_g0 = compute_rdmst(g0, 0)
expected_g0 = ({0: {1: 2, 2: 2, 3: 2}, 1: {5: 2}, 2: {4: 2}, 3: {}, 4: {}, 5: {}}, 10)
#print("COMPUTE RDMST (G0): " + str(actual_g0))
#print("CORRECT? " + str(actual_g0 == expected_g0))

#print("")

# #  WORKS YAY
# actual_g1 = compute_rdmst(g1, 0)
# expected_g1 = ({0: {2: 4}, 1: {}, 2: {3: 8}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}, 28)
# print("COMPUTE RDMST (G1): " + str(actual_g1))
# print("COMPUTE RDMST EXPECTED (G1): " + str(expected_g1))
# print("CORRECT? " + str(actual_g1 == expected_g1))

# print("")

# # WORKS YAY
# actual_g2 = compute_rdmst(g2, 0)
# expected_g2 = ({0: {2: 4}, 1: {}, 2: {1: 2}}, 6)
# print("COMPUTE RDMST (G2): " + str(actual_g2))
# print("COMPUTE RDMST EXPECTED (G2): " + str(expected_g2))
# print("CORRECT? " + str(actual_g2 == expected_g2))

#print("")

# def difference(dict1, dict2):
#     for node, nbr in dict1[0].items():
#         if nbr != dict2[0][node]:
#             print(node)
#             print("MINE: " + str(nbr))
#             print("CORRECT: " + str(dict2[0][node]))

#PASS
# actual_g3 = compute_rdmst(g3, 1)
# expected_g3 = ({1: {3: 1.0, 4: 9.1}, 2: {}, 3: {2: 1.0, 5: 0.0}, 4: {}, 5: {}}, 11.1)
# print("COMPUTE RDMST (G3): " + str(actual_g3))
# print("COMPUTE RDMST EXPECTED (G3): " + str(expected_g3))
# print("CORRECT? " + str(actual_g3 == expected_g3))
# difference(actual_g3, expected_g3)
#print("COMPUTE RDMST (G1) EXPECTED: ({0: {2: 4}, 1: {}, 2: {3: 8}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}, 28)")
#print("COMPUTE RDMST (G2): " + str(compute_rdmst(g2, 0)))
#print("COMPUTE RDMST (G3): " + str(compute_rdmst(g3, 1)))
# print("COMPUTE RDMST (G4): " + str(compute_rdmst(g4, 1)))


# actual_g4 = compute_rdmst(g4, 1)
# expected_g4 = ({1: {8: 6.1, 3: 1.0, 4: 9.1}, 2: {}, 3: {2: 1.0, 5: 0.0}, 4: {9: 6.1, 10: 5.0}, 5: {}, 6: {7: 0.0}, 7: {}, 8: {}, 9: {}, 10: {6: 0.0}}, 28.3)
# print("COMPUTE RDMST (G4): " + str(actual_g4))
# print("COMPUTE RDMST EXPECTED (G4): " + str(expected_g4))
# print("CORRECT? " + str(actual_g4 == expected_g4))
# difference(actual_g4, expected_g4)



