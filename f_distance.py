#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:41:42 2019

@author: Dimos Goundaroulis (Gkountaroulis)

This program computes the experimental f-distance of isotopy classes of 
knotoid diagrams, as presented in the paper:

f-distance of knotoids and protein structure by A. Barbensi and D. Goundaroulis
(arxiv preprint 2019)

The program uses the classification data for knotoids on the sphere from:
A systematic classification of knotoids on the plane and on the sphere by:
D. Goundaroulis, J. Dorier, M. Lackenby and A. Stasiak (Submitted 2019. See also arXiv:1902.07277).

To compute the experimental f-distance, we create a graph where the vertices correspond to 
knotoid diagrams. An edge between two vertices means that the corresponding diagrams are related
with a single forbidden move. 

The experimental f-distances are computed by partitioning the set of all diagrams into isotopy classes and
by finding the Dijkstra path between all pairs of classes.
"""

import pickle
import pandas
from optparse import OptionParser
from collections import deque
import igraph


'''

Below are two auxilliary functions for 
searching in lists and in lists of lists, respectively.

'''
def find_in_list(mylist,target):
    for i in range(len(mylist)):
        if mylist[i]==target:
            yield i

def find_in_list_of_lists(mylist_of_lists,target):
    for i in range(len(mylist_of_lists)):
        if target in mylist_of_lists[i]:
            yield i



def str_to_gcode(gcodestr):
    '''

    Converting a Gauss code in string format to a list. 
    The function detects if the Gauss code corresponds to spherical or planar knotoid
    and it treats it accordingly.

    '''
    gcode = gcodestr.split('\t')
    gcode[0] = [int(x) for x in gcode[0].split(' ')]
    help = []
    for i in range(len(gcode[1])):
        if gcode[1][i] == '+':
            help.append(1)
        if gcode[1][i] == '-':
            help.append(-1)
    gcode[1] = help
    if len(gcode)==3:
        gcode[2] = [int(x) for x in gcode[2].split(' ')]
    return gcode



def create_gc(gc):
    '''

    Function for creating a Gauss Code

    '''
    global flag_planar
    if type(gc)==str:
        gcode = str_to_gcode(gc)
    if gcode:
        try:
            flag_planar=True
            return [gcode[0],gcode[1],gcode[2]]
        except IndexError:
            flag_planar=False
            return [gcode[0],gcode[1]]


def convert_to_pd(gc):
    '''

    Convert Gauss Code to PD code


    '''
    PD = []
    checked = []
    for x in range(len(gc[0])):
        inv = next(find_in_list(gc[0],-1 * (gc[0][x])))
        if abs(gc[0][x]) not in checked:
            if gc[0][x] > 0 and gc[1][abs(gc[0][x]) - 1] == 1:
                PD.append([inv, x + 1, inv + 1, x])
            if gc[0][x] > 0 and gc[1][abs(gc[0][x]) - 1] == -1:
                PD.append([inv, x, inv + 1, x + 1])
            if gc[0][x] < 0 and gc[1][abs(gc[0][x]) - 1] == 1:
                PD.append([x, inv + 1, x + 1, inv])
            if gc[0][x] < 0 and gc[1][abs(gc[0][x]) - 1] == -1:
                PD.append([x, inv, x + 1, inv + 1])
            checked.append(abs(gc[0][x]))
    return PD



def regions(inpt):
    '''

    Computes the local regions of a knotoid diagram. A region is bounded by arcs of the diagram.
    The function starts from the first crossing and moves in a counter-clockswise fashion following
    the arcs that bound the region, until it returns to the starting point. Each arc is assigned a
    sign: +1 if the orientation of the search agrees with the orientation of the arc and -1 if it 
    doesn't agree. Note that the search algorithm goes around an endpoint and so the adjacent arc will
    appear twice (but with different signs) in the corresponding region.

    '''
    checked = set()
    region = set()
    regionlist = []
    nb_endpoints = 0
    pd=convert_to_pd(inpt)
    for i in range(len(pd)):
        for k in range(4):
            if (i, k) not in checked:
                region = []
                i1 = i
                k1 = k
                while True:
                    checked.add((i1, k1))
                    if pd[i1][k1] - pd[i1][(k1 + 2) % 4]==1 or pd[i1][k1] - pd[i1][(k1 + 2) % 4]<-1:
                        region.append((pd[i1][k1], 1))
                    elif pd[i1][k1] - pd[i1][(k1 + 2) % 4]==-1 or pd[i1][k1] - pd[i1][(k1 + 2) % 4]>1:
                        region.append((pd[i1][k1], -1))

                    found = False
                    for j in range(len(pd)):
                        for l in range(4):
                            if pd[i1][k1] == pd[j][l] and not (i1 == j and k1 == l):
                                i1 = j
                                k1 = (l + 1) % 4
                                found = True
                                break
                        if found:
                            break
                    if not found:
                        nb_endpoints += 1
                        if pd[i1][k1] > pd[i1][(k1 + 2) % 4]:
                            region.append((pd[i1][k1], -1))
                        else:
                            region.append((pd[i1][k1], 1))
                        if nb_endpoints > 2:
                            print("too many endpoints")
                            return []
                        k1 = (k1 + 1) % 4
                    if i1 == i and k1 == k:
                        break
                regionlist.append(region)
    return regionlist



def sign_outside(gc, the_regions):
    '''
    
    Detection of the outer (unbounded) region of the diagram.
    
    '''
    outside_region = []
    for region in the_regions:
        if set(gc[2]) == set([m[0] for m in region]):
            outside_region = region

    return outside_region


def renumber(temp):
    '''
    
    Function for renumbering of a Gauss code.
    
    
    '''
    checked = []
    for i in range(0, len(temp)):
        if abs(temp[i]) not in checked:  
            checked.append(abs(temp[i]))
            if temp[i] < 0:  
                temp[i] = -(checked.index(abs(temp[i])) + 1)
            else:
                temp[i] = checked.index(abs(temp[i])) + 1
        else:
            if temp[i] < 0:
                temp[i] = -(checked.index(abs(temp[i])) + 1)
            else:
                temp[i] = checked.index(abs(temp[i])) + 1
    return temp


def tail_forbidden(gc):
    '''
    
    Performs a forbidden move at the TAIL of a knotoid diagram. This function effectively removes
    a crossing from the corresponding Gauss code and calls the renumbering function on the resulting code.
    
    '''
    newgc=pickle.loads(pickle.dumps(gc, -1))

    if flag_planar:
        pd = convert_to_pd(newgc)

        crossing_with_endpoint = pd[next(find_in_list_of_lists(pd,0))]
        min_arc_idx = next(find_in_list(crossing_with_endpoint,0))
        other_strand=sorted((crossing_with_endpoint[(min_arc_idx-1)%4], crossing_with_endpoint[(min_arc_idx+1)%4]))
        first=next(find_in_list(newgc[0],1))
        newgc[0].pop(first)
        second=next(find_in_list(newgc[0],-1))
        newgc[0].pop(second)
        newgc[1].pop(0)
        newgc[0]=renumber(newgc[0])

        if newgc[0]==[] and newgc[1]==[]:
            newgc[2]=[0]
            
        elif (newgc[0]==[] and newgc[1]!=[]) or (newgc[0]!=[] and newgc[1]==[]):
            raise Exception("Error in Gauss Code.")
            

        
        elif 0 in newgc[2] and all(x in newgc[2] for x in other_strand):

            if 1 in newgc[2]:
                the_regions=regions(newgc)
                for region in the_regions:

                    if (0,1) in region and (0,-1) in region:

                        newgc[2]=[x[0] for x in region]
                        break
                
            else:
                for i in range(len(newgc[2])):

                    if newgc[2][i]<=other_strand[0]:
                        newgc[2][i]=abs(newgc[2][i]-1)
                    elif newgc[2][i]>=other_strand[1]:
                        newgc[2][i]=abs(newgc[2][i]-2)
                newgc[2].pop(0)

            
        
        
        elif 0 not in newgc[2] and any(x in newgc[2] for x in other_strand):

            the_regions=regions(newgc)

            if the_regions:
                for region in the_regions:

                    if (0,1) in region and (0,-1) in region:
                        newgc[2]=[x[0] for x in region]
                        break
            else:
                newgc[2]=[0]

        elif 0 not in newgc[2] and not any(x in newgc[2] for x in other_strand):
            for i in range(len(newgc[2])):
                if newgc[2][i]<=other_strand[0]:
                    newgc[2][i]=newgc[2][i]-1
                elif newgc[2][i]>=other_strand[1]:
                    newgc[2][i]=newgc[2][i]-2

    
        newgc[2]=sorted(list(set(newgc[2])))
    else:
        first=next(find_in_list(newgc[0],1))
        newgc[0].pop(first)
        second=next(find_in_list(newgc[0],-1))
        newgc[0].pop(second)
        newgc[1].pop(0)
        newgc[0]=renumber(newgc[0])

    return newgc
            
        
    

def head_forbidden(gc):
    '''
    
    Performs a forbidden move at the HEAD of a knotoid diagram. This function effectively removes
    a crossing from the corresponding Gauss code and calls the renumbering function on the resulting code.
    
    '''
    newgc=pickle.loads(pickle.dumps(gc, -1))
    if flag_planar:
        pd = convert_to_pd(newgc)
        max_arc=2*max(newgc[0])
        crossing_with_endpoint = pd[next(find_in_list_of_lists(pd,max_arc))]
        max_arc_idx = next(find_in_list(crossing_with_endpoint,max_arc))
        other_strand=(crossing_with_endpoint[(max_arc_idx-1)%4], crossing_with_endpoint[(max_arc_idx+1)%4])

        first=newgc[0][-1]
        newgc[0].pop()
        second=next(find_in_list(newgc[0],-first))
        newgc[0].pop(second)
        newgc[0]=renumber(newgc[0])
        newgc[1].pop(abs(first)-1)
        
        if newgc[0]==[] and newgc[1]==[]:
            newgc[2]=[0]
            
        elif (newgc[0]==[] and newgc[1]!=[]) or (newgc[0]!=[] and newgc[1]==[]):
            raise Exception("Error in Gauss Code.")
            
        elif max_arc in newgc[2] and all(x in newgc[2] for x in other_strand):
            if max_arc-1 in newgc[2]:
                the_regions=regions(newgc)
                for region in the_regions:
                    if (max_arc-2,1) in region and (max_arc-2,-1) in region:
                        newgc[2]=[x[0] for x in region]
                        break
                
            else:
                largest_arc = max(other_strand)
                newgc[2] = [x-1 if x>=largest_arc else x for x in newgc[2]]
                newgc[2].pop()


        elif max_arc not in newgc[2] and any(x in newgc[2] for x in other_strand):
            the_regions=regions(newgc)
            if the_regions:
                for region in the_regions:
                    if (max_arc-2,1) in region and (max_arc-2,-1) in region:
                        newgc[2]=[x[0] for x in region]
                        break
            else:
                newgc[2]=[0]

        elif max_arc not in newgc[2] and not any(x in newgc[2] for x in other_strand):
            largest_arc = max(crossing_with_endpoint[(max_arc_idx-1)%4], crossing_with_endpoint[(max_arc_idx+1)%4])
            for i in range(len(newgc[2])):
                if newgc[2][i] >= largest_arc:
                    newgc[2][i] = newgc[2][i]-1
        newgc[2]=sorted(list(set(newgc[2])))
    else:
        first=newgc[0][-1]
        newgc[0].pop()
        second=next(find_in_list(newgc[0],-first))
        newgc[0].pop(second)
        newgc[0]=renumber(newgc[0])
        newgc[1].pop(abs(first)-1)

    return newgc



def add_edge_data_forbidden(g):
    '''
    
    Adds the edges in the graph of knotoid diagrams. 
    An edge connects two diagrams that are related with a single forbidden move.
    
    '''
    tail_f = tail_forbidden(g)
    edges_to_add.append((vertex_to_index_list[repr(g)],vertex_to_index_list[repr(tail_f)]))

    head_f = head_forbidden(g)
    edges_to_add.append((vertex_to_index_list[repr(g)],vertex_to_index_list[repr(head_f)]))

  
    

'''
    Main
'''
  
     
if __name__ == "__main__":

    parser = OptionParser()   
    parser.add_option("-f", "--file", help='Specify input filename.', action="store",type="string",dest="filename")
    parser.add_option("-c", "--composite", help="Consider composites in the analysis.", action="store_true", default=False, dest="composites")
    parser.add_option("-p", "--plot", help="Plots connectivity graph.", action = "store_true", default=False, dest = "plotting")
    parser.add_option("--debug", help='Debug mode.', action='store_true', default=False, dest='debug')
    parser.add_option("--start_part", help="Partition tracking. Starting partition. Use only with --debug mode.", action="store", dest="start_part")
    parser.add_option("--end_part", help="Partition tracking. Ending partition. Use only with --debug mode.", action="store", dest="end_part")
    (options,args)=parser.parse_args()
    optionsdict = vars(options)
    inputfilename=options.filename
    composites=options.composites
    plotting=options.plotting
    debug=options.debug
    start_part=options.start_part
    end_part=options.end_part

    if inputfilename==None:
        raise Exception('No input file specified.')
        

    print('Performing forbidden moves.')
        
    if composites:
        print ('Including composites.')
    else:
        print ('Not including composites.')
    
    if debug:
        if (start_part and not end_part) or (not start_part and not end_part) or (not start_part and not end_part):
            raise Exception ("Debug must be used with bot --start_part and --end_part.")

    if not debug:
        if start_part or end_part:
            raise Exception ("Options --start_part and --end_part must be used together with --debug.")


    file = open(options.filename, 'r')
    gdict=deque()
    vertex_to_index_list=dict()

    isotopy_vertices=set()
    check=file.readline().strip('\n').split('\t')
    me=1
    for l in file:
        line = l.strip('\n').split('\t')
        if check[2]=="gc3":
            gc = create_gc('\t'.join(line[0:3]))
            status=line[3]
            iso_class=line[4]
            gdict.append([gc,status,iso_class])
            if not composites:
                if status!="composite":
                    isotopy_vertices.add(iso_class)
            else:
                isotopy_vertices.add(iso_class)
            vertex_to_index_list[repr(gc)]=me
            me+=1
        else:
            gc = create_gc('\t'.join(line[0:2]))
            status=line[2]
            iso_class=line[3]
            gdict.append([gc,status,iso_class])
            if not composites:
                if status!="composite":
                    isotopy_vertices.add(iso_class)
            else:
                isotopy_vertices.add(iso_class)
            vertex_to_index_list[repr(gc)]=me
            me+=1
            
            
    file.close()
    tmp=deque(x[0] for x in gdict)
    if flag_planar:
        tmp.appendleft([[],[],[0]])
        vertex_to_index_list[repr([[],[],[0]])]=0
    else:
        tmp.appendleft([[],[]])
        vertex_to_index_list[repr([[],[]])]=0


    print("Creating graph.")

    graph=igraph.Graph(n=len(gdict)+1,directed=False)
    
    if flag_planar:
        gdict.appendleft([[[[],[],[0]]],"prime","0_1"])
    else:
        gdict.appendleft([[[[],[]]],"prime","0_1"])
    
    print("Adding vertex attributes.")
    graph.vs["gc"] =[x[0] for x in gdict]
    graph.vs["status"] = [x[1] for x in gdict]
    graph.vs["iso_class"] = [x[2] for x in gdict]
    

        
    
    print("Creating edge data.")

    file_name=options.filename.strip('.txt')
     
    edges_to_add=deque()
    for el in tmp:
        if flag_planar:
            if el==[[],[],[0]]:
                continue
        else:
            if el==[[],[]]:
                continue
        add_edge_data_forbidden(el)


    print("Adding edges to graph.")
    try:
        graph.add_edges(edges_to_add)
    except NameError:
        raise Exception("Vertex doesn't exist.")

    if debug:
        print("Moves from {} to {}".format(start_part, end_part))
        flname="Moves_from_{}_to_{}.txt".format(start_part,end_part)
        file=open(flname,"w")
        file.write("From isotopy class"+"\t"+"Diagram"+"\t"+"To isotopy class"+"\t"+"Diagram"+"\n")
           
        for e in graph.es:
            if start_part in graph.vs[e.source]["iso_class"] and end_part in graph.vs[e.target]["iso_class"]:
                file.write(str(graph.vs[e.source]["iso_class"])+"\t"+str(graph.vs[e.source]["gc"])+"\t"+str(graph.vs[e.target]["iso_class"])+"\t"+str(graph.vs[e.target]["gc"])+"\n")
        file.close()
                

    print("-------------------")
    print("Full graph summary.")       
    igraph.summary(graph)  
    connectivity_list=[]
    print("-------------------")
#
    if not composites:
        to_delete_ids = [v.index for v in graph.vs if 'ICID' in v['iso_class']]
        graph.delete_vertices(to_delete_ids)
    
        print("Creating isotopy graph.")
        for e in graph.es:
            if 'ICID' not in graph.vs[e.source]["iso_class"] or 'ICID' not in graph.vs[e.target]["iso_class"]:
                connectivity_list.append((graph.vs[e.source]["iso_class"],graph.vs[e.target]["iso_class"]))
    else:
        for e in graph.es:
            connectivity_list.append((graph.vs[e.source]["iso_class"],graph.vs[e.target]["iso_class"]))
    
    
    isotopy_graph=igraph.Graph(n=len(isotopy_vertices),directed=False)
    
    j=0
    new_vertex_name_to_index_list=dict()

    for v in sorted(list(isotopy_vertices)):
        new_vertex_name_to_index_list[v]=j
        isotopy_graph.vs[j]["label"]=v
        j+=1


    new_edge_list=set()
    for pair in connectivity_list:
        new_edge_list.add((new_vertex_name_to_index_list[pair[0]],new_vertex_name_to_index_list[pair[1]]))
    
    isotopy_graph.add_edges(new_edge_list)
    simplified_graph = isotopy_graph.simplify(multiple=True, loops=True)
    

    print("Creating distance matrix.")
    dist_matrix=simplified_graph.shortest_paths_dijkstra(source=simplified_graph.vs, target=simplified_graph.vs, weights=None)
    df = pandas.DataFrame(dist_matrix,index=simplified_graph.vs['label'], columns = simplified_graph.vs['label'])
    df = df.sort_index()
    df = df.reindex(sorted(df.columns), axis=1)
    df.to_csv(file_name+"_forbidden_data.txt",index=True,sep="\t")


    if plotting:
        print("Plotting graph.")
        layout = simplified_graph.layout("kk")
        visual_style = {}
        visual_style["vertex_size"] = 40
        visual_style["layout"] = layout
        visual_style["bbox"] = (3840, 2160) 
        visual_style["margin"] = 100
        out=igraph.plot(simplified_graph, **visual_style)
        out.save(file_name+"_forbidden_graph.png")

    print("Done.")
