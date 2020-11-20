#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 10:15:26 2019

@author: dgkounta
"""
from optparse import OptionParser
from subprocess import Popen, DEVNULL,STDOUT
import igraph
import pandas
import shlex
import multiprocessing as mp
from scipy.spatial import SphericalVoronoi
import numpy
from spherical_geometry.polygon import SphericalPolygon


'''

   Main Function
   Protein Analysis
   
'''
def load_pdb_file(name):
    file=open(name+".pdb","r")
    protdata=file.readlines()
    file.close()
    return protdata

def get_backbone_coordinates(protdata,chain):
    backbone=[]
    for i in range(len(self.protdata)):
            line=str.split(self.protdata[i])
            # print (line[0])
            if(line[0]=="ATOM" and line[2]=="CA" and self.protdata[i][21]==chain):

                x=float(self.protdata[i][30:38])
                y=float(self.protdata[i][38:46])
                z=float(self.protdata[i][46:54])
                resid=int(self.protdata[i][22:26])
                backbone.append([x,y,z,resid])

            if line[0]=="ENDMDL":
                break
    # print (backbone)
    return np.array(backbone)



def analyze_protein(protein):
    prot=protein[0]
    chain=protein[1]
    print ("Downloading {}{}.".format(prot,chain))

    cmd="wget https://files.rcsb.org/download/{}.pdb".format(prot)
    proc=Popen(shlex.split(cmd), stdout=DEVNULL, stderr=STDOUT)
    proc.communicate()
    

    print ("Extracting coordinates.")
    cmd="./pdb_to_xyz.R {}.pdb --chain-id={} -o {}.txt".format(prot,chain,prot)
    proc=Popen(shlex.split(cmd), stdout=DEVNULL, stderr=STDOUT)
    proc.communicate()


    print ("Removing pdb file.")
    proc=Popen(["rm","{}.pdb".format(prot)],stdout=DEVNULL, stderr=STDOUT)
    proc.communicate()

    print ("Analyzing protein")





    for proj in projections:
        cmd="./polynomial_invariant {}.txt --names-db=internal --arrow --nb-projections={} --output-diagram={}{}_gl_{}.txt".format(prot,proj,prot,chain,proj)
        proc=Popen(shlex.split(cmd), stdout=DEVNULL, stderr=STDOUT)
        proc.communicate()
    
        name="{}{}_gl_{}.txt".format(prot,chain,proj)
        origins,origin_types=[],[]
        file=open(name,"r")        
        next(file)
        me=0
        checked=[]
        vertex_to_index_list=dict()
        edges_to_add=[]
        for l in file:
            line = l.strip('\n').split('\t')
            if line[3]!='UNKNOWN':
                if '*' not in line[3]:
                    origin_types.append(line[3])
                    origins.append([float(line[0]),float(line[1]),float(line[2])])
                    if line[3] not in checked:
                        vertex_to_index_list[repr(line[3])]=me
                        me+=1
                
            
        origins=numpy.array(origins)
        voro=SphericalVoronoi(origins)
        voro.sort_vertices_of_regions()
        for i in range(len(voro.regions)):
            for j in range(i+1 ,len(voro.regions)):
                intersect=list(set(voro.regions[i]) & set(voro.regions[j]))
                if intersect:
                    permissable = set([1,len(voro.regions[i]),len(voro.regions[j])])
                    index_A = [voro.regions[i].index(x) for x in intersect]
                    diff_A = [abs(j-i) for i, j in zip(index_A[:-1], index_A[1:])]
                    index_B = [voro.regions[j].index(x) for x in intersect]
                    diff_B = [abs(j-i) for i, j in zip(index_B[:-1], index_B[1:])]
                    # print ("intersect",intersect,"voro[i]",voro.regions[i],"voro[j]",voro.regions[j])
                    # print ("index_A",index_A,"diff_A",diff_A,"index_B",index_B,"diff_B",diff_B)
                    # print ("perm_A",permissable & set(diff_A),"perm_B", permissable & set(diff_B))
                    if permissable & set(diff_A) and permissable & set(diff_B):
                        for k in range(int(len(intersect)/2)):
                            edges_to_add.append((vertex_to_index_list[repr(origin_types[i])],vertex_to_index_list[repr(origin_types[j])]))

# Comment lines 155 - 165 and uncomment lines 169-170 to go to previous version.

                    # for k in range(int(len(intersect)/2)):
                    #     edges_to_add.append((vertex_to_index_list[repr(origin_types[i])],vertex_to_index_list[repr(origin_types[j])]))
      
        # print ("len",len(voro.vertices))
        print("Creating graph.")
    
        g=igraph.Graph(len(origin_types),directed=False)
        
        
        origin_types=[x.split("|")[0].strip('m').strip('s') for x in origin_types]
        g.vs["label"]=[x for x in origin_types]
        g.add_edges(edges_to_add)
        all_connections=[]
        for edge in g.es:
          src = g.vs[edge.source]['label']
          tgt = g.vs[edge.target]['label']
          if (src=='3_1' and tgt!='3_1'):
              all_connections.append((src,tgt))
          elif  (src!='3_1' and tgt=='3_1'):
              all_connections.append((tgt,src))

        
        proj_errors=0
        error_pairs=[]
        for pair in all_connections:

            if theoretical_values[pair[0]][pair[1]]!=1:
                proj_errors+=1
                error_pairs.append(pair)
        try:
            errors=round(proj_errors/len(all_connections),4)
        except ZeroDivisionError:
            errors='NaN'
            


        edges_between=[]
        for edge in g.es:
          src = g.vs[edge.source]['label']
          tgt = g.vs[edge.target]['label']
          if (src=='3_1' and tgt=='2_1'):
              edges_between.append((src,tgt))
          elif  (src=='2_1' and tgt=='3_1'):
              edges_between.append((tgt,src))

        if len(edges_between)>0:
            try:
                interface=round(len(edges_between)/len(all_connections),4)
            except ZeroDivisionError:
                interface='NaN'
        else:
            interface='NaN'
        
        area_31,area_21=[],[]
        all_area_31,all_area_21=[],[]

        # voro_areas = voro.calculate_areas()
        
        print ("voro_areas: ",voro_areas)
        for i in range(len(voro.regions)):
            voronoi_cell=numpy.asarray([list(voro.vertices[x]) for x in voro.regions[i]])
            if origin_types[i]=='3_1':
                sph_pol=SphericalPolygon(voronoi_cell)
                area_31.append(sph_pol.area())
            elif origin_types[i]=='2_1':
                sph_pol=SphericalPolygon(voronoi_cell)

                area_21.append(sph_pol.area())
        all_area_31=numpy.sum(area_31)
        all_area_21=numpy.sum(area_21)
        

        all_graphs.append({'n':str(proj),'Error':errors, 'Interface': interface, 'Protein':prot+chain, 'Spectrum':len(set(origin_types)),
                           'Error Pairs':set(error_pairs),'Full Spectrum':set(origin_types),
                           'Area  3_1':all_area_31, 'Area 2_1':all_area_21,'Area 2_1/3_1': all_area_21/all_area_31})
    
    
if __name__ == "__main__":

    parser = OptionParser()   
    parser.add_option("-f", "--file", help='Specify input filename with protein names.', action="store",type="string",dest="filename")
    parser.add_option("-n", "--nb-projections", help='Number of projections.', action="store", type="string", default=100, dest='projections')
    parser.add_option("-p", "--plot", help='Plot graphs.', action='store_true', default=False, dest='plot')
    parser.add_option("-i", "--input-protein", help='PDB code to analyze.', action="store", type="string",dest='inputprotein')
    parser.add_option("-c", "--chain", help='Protein chain to analyze.', action='store', type='string',default='A',dest='chain')

    (options,args)=parser.parse_args()
    optionsdict = vars(options)
    inputfilename=options.filename
    projections=options.projections
    plotting=options.plot
    inputprotein=options.inputprotein
    chain=options.chain

    try:
        projections=[int(x) for x in projections.split(",")]
    except AttributeError:
        projections=[projections]

    if inputfilename==None and inputprotein==None:
        raise Exception('No input file or input PDB code specified.')    
        
    if inputfilename and inputprotein:
        raise Exception("The -i option can't be used together with the -f option.")


    print ('Loading theoretical values.')
    global theoretical_values
    theoretical_values=pandas.read_csv('theoretical_values.csv',sep="\t",index_col=0)
    
    if inputfilename:
        with open(inputfilename,"r") as f:
            file=f.readlines()
        
        proteins=[]
        manager=mp.Manager()
    
        all_graphs=manager.list()
        for l in file[1:]:
            line=l.strip("\n").strip('  ').split("|")  ### Old version
            # line = l.strip("\n").split(";") ### knotprot api file
            print ("lll",line)
            # if line[4]=='K' and line[5]=='31':
            proteins.append((line[0],line[1]))
        # file.close()
        pool = mp.Pool(mp.cpu_count())
        for i in range(len(proteins)):
#            analyze_protein(proteins[i])
            pool.apply_async(analyze_protein,(proteins[i],))
        pool.close()
        pool.join()
        all_graphs=list(all_graphs)

    if inputprotein:
        all_graphs=[]
        analyze_protein([inputprotein,chain])



    print('Creating Dataframe.')
    df = pandas.DataFrame(all_graphs)
    if plotting:
        df=df.set_index(['n'])
        ax = df.plot(kind='line',subplots=True,title=inputprotein)[0].get_figure()
        if inputprotein:
            ax.savefig('graph_{}.pdf'.format(inputprotein))
        else:
            ax.savefig('graph_all_proteins.pdf')

    print('Saving Dataframe.')
    if plotting:
        df=df.reset_index()
    print ("df:", df)
    df=df.set_index(['Protein','n']).sort_index(level=0)
    if inputprotein:
        df.to_csv('Protein_data_{}.csv'.format(inputprotein),sep="\t")
    else:
        df.to_csv('Protein_data_all.csv',sep="\t")

    

    print ('Done.')
        

    
    
    
            
# 