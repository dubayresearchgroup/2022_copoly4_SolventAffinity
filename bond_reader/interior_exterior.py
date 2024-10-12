import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import chainanalysis
import os
import sys
import subprocess as sb
import pickle
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors




def calculate_int_ext(filename,criteria_int,criteria_ext):
    #Filenames and properties for I/O
    #filename = '200813_copoly3_standard-1'
    #criteria_int = 60
    #criteria_ext = 12
    Nfilter = 5 #minimum number of monomers in cluster for splitting
    inputfile = open(filename+'.data')
    inputlines = inputfile.readlines()
    numheaderlines = 33 #DONT TOUCH THIS: 33 for standard 21 for early output
    outputfile = open('VMD.data','w')
    tempfile = open('temp.csv','w')
    typefile = open(filename+'-type')
    typelines = typefile.readlines()
    typefile.close()
    numatoms = int(str.split(inputlines[2])[0])
    numbonds = int(str.split(inputlines[4])[0])
    numangles = int(str.split(inputlines[6])[0])
    numskips = 3



    #Splitting out atom coordinates for clustering and writing correct data format snapshot for VMD
    try:
        tempfile.write('atomid,type,x,y,z\n')
        xhi = float(str.split(inputlines[9])[1]) #for unwrapping pbc coordinates (doesn't work well)
        yhi = float(str.split(inputlines[10])[1])
        zhi = float(str.split(inputlines[11])[1])
        numatoms = int(str.split(inputlines[2])[0])
        numbonds = int(str.split(inputlines[4])[0])
        numangles = int(str.split(inputlines[6])[0])
        numskips = 3
        for line in inputlines[:numheaderlines]:
            outputfile.write(line)
        for line in inputlines[numheaderlines:]:
            split = str.split(line)
            if len(split)>7 and len(split)<11:
                outputfile.write(split[0]+' '+split[2]+' '+split[2]+' '+split[3]+' '+split[4]+' '+split[5]+'\n')
                tempfile.write(split[0]+','+split[2]+','+split[3]+','+split[4]+','+split[5]+'\n')
            else:
                outputfile.write(line)
            
    finally:
        inputfile.close()
        outputfile.close()
        tempfile.close()
        
     
     
    #Clustering analysis 
    print('Generating clusters.')
    df = pd.read_csv('temp.csv')
    df = df.sort_values(by=['atomid'])
    df_centers = df.loc[df['type']>1]
    X_centers = df_centers.values
    neighbors_model = NearestNeighbors(radius = 2.5)
    neighbors_model.fit(X_centers[:,2:])
    neighborhoods = neighbors_model.radius_neighbors(X_centers[:,2:],return_distance=False)
    X_int = []
    int_centers = []
    int_ids = []
    X_ext = []
    ext_centers = []
    ext_ids = []
    X_both = []
    both_centers = []
    both_ids = []
    num_neighbors = []
    for i,neighbor_list in enumerate(neighborhoods):
        num_neighbors.append(len(neighbor_list))
        if len(neighbor_list) > criteria_int:
            X_int.append(X_centers[i,:].tolist())
            int_centers.append(int(X_centers[i,0]))
            X_both.append(X_centers[i,:].tolist())
            both_centers.append(int(X_centers[i,0]))
        elif len(neighbor_list) > criteria_ext and criteria_ext <= 60: 
            X_ext.append(X_centers[i,:].tolist())
            ext_centers.append(int(X_centers[i,0]))
            X_both.append(X_centers[i,:].tolist())
            both_centers.append(int(X_centers[i,0]))
            
    X_int = np.array(X_int)
    X_ext = np.array(X_ext)
    X_both = np.array(X_both)

    mono_int = len(X_int[:,1])
    mono_ext = len(X_ext[:,1])
    mono_both = len(X_both[:,1])
    int_As = 0
    int_Bs = 0
    ext_As = 0
    ext_Bs = 0
    both_As = 0
    both_Bs = 0
    for center in X_int[:,1]:
        if int(center)==2:
            int_As += 1
        else:
            int_Bs += 1
            
    for center in X_ext[:,1]:
        if int(center)==2:
            ext_As += 1
        else:
            ext_Bs += 1
            
    for center in X_both[:,1]:
        if int(center)==2:
            both_As += 1
        else:
            both_Bs += 1

    FA_int = int_As/mono_int
    FB_int = int_Bs/mono_int
    FA_ext = ext_As/mono_ext
    FB_ext = ext_Bs/mono_ext
    FA_both = both_As/mono_both
    FB_both = both_Bs/mono_both

    print('Internal FA {}'.format(FA_int))
    print('Internal FB {}'.format(FB_int))
    print('External FA {}'.format(FA_ext))
    print('External FB {}'.format(FB_ext))
    print('Full structure FA {}'.format(FA_both))
    print('Full structure FB {}'.format(FB_both))

    monos_ratio = [int_As,int_Bs,mono_int,ext_As,ext_Bs,mono_ext,both_As,both_Bs,mono_both]

    with open('./jar/'+filename+'_int-ext.pickle','wb') as f:
        pickle.dump(monos_ratio,f)


    for center in int_centers:
        int_ids.append(center-1)
        int_ids.append(center)
        int_ids.append(center+1)
    for center in ext_centers:
        ext_ids.append(center-1)
        ext_ids.append(center)
        ext_ids.append(center+1)
    for center in both_centers:
        both_ids.append(center-1)
        both_ids.append(center)
        both_ids.append(center+1)
        
    os.remove('temp.csv')



    #plotting
    fig= plt.figure(figsize=(16,16),dpi=300)
    ax = fig.add_subplot(221, projection='3d')


    ax1 = fig.add_subplot(222, projection='3d')
    ax2 = fig.add_subplot(223, projection='3d')
    ax3 = fig.add_subplot(224, projection='3d')

    ax.scatter(X_centers[:,2], X_centers[:,3],X_centers[:,4],c='k',s=0.5)
    ax1.scatter(X_both[:,2], X_both[:,3],X_both[:,4],c='b',s=0.5)
    ax2.scatter(X_int[:,2], X_int[:,3],X_int[:,4],c='r',s=0.5)
    ax3.scatter(X_both[:,2], X_both[:,3],X_both[:,4],c='b',s=0.5,alpha=0.25)
    ax3.scatter(X_int[:,2], X_int[:,3],X_int[:,4],c='r',s=0.5,alpha=1.)
    plt.savefig(filename+'_interior_exterior_plot.png')








    #filtering
    for val in [0,1,2]:
        if val == 0:
            ids = int_ids
            id_strings = [str(id) for id in int_ids]
            string = '_internal'+'_cut'+str(criteria_int)
        elif val == 1:
            ids = ext_ids
            id_strings = [str(id) for id in ext_ids]
            string = '_external'+'_cut'+str(criteria_ext)
        else:
            ids = both_ids
            id_strings = [str(id) for id in both_ids]
            string = '_both'



        outputfile = open(filename+string+'.data','w')
        clusteratoms = 0
        clusterbonds = 0
        clusterangles = 0
        entrynum = 0
        entries = []
        print('Splitting '+string+' files.')


        #THIS BLOCK SHOULD SKIP SMALL CLUSTERS - ADJUST SIZE BY NFILTER
        if len(id_strings) < (3*Nfilter):
            print('Cluster less than {} monomers, ignoring.'.format(Nfilter))
            outputfile.close()

            #Cleanup
            os.remove(filename+string+'.data')


        else:

                #Data file
            print('Splitting data file.')


            #Header
            for line in inputlines[:numheaderlines]:
                outputfile.write(line)

            #Atoms
            for line in inputlines[numheaderlines:numheaderlines+numatoms]:
                split = str.split(line)
                if split[0] in id_strings:
                    outputfile.write(split[0]+' '+split[2]+' '+split[2]+' '+split[3]+' '+split[4]+' '+split[5]+'\n')
                    clusteratoms += 1
            #Velocities
            for line in inputlines[numheaderlines+numatoms:numheaderlines+numatoms+numskips]:
                outputfile.write(line)
            for line in inputlines[numheaderlines+numatoms+numskips:2*numatoms+numheaderlines+numskips]:
                split = str.split(line)
                if split[0] in id_strings:
                    outputfile.write(line)

            #Bonds
            for line in inputlines[2*numatoms+numheaderlines+numskips:2*numatoms+numheaderlines+2*numskips]:
                outputfile.write(line)
            for line in inputlines[2*numatoms+numheaderlines+2*numskips:2*numatoms+numheaderlines+2*numskips+numbonds]:
                split = str.split(line)
                if split[2] in id_strings and split[3] in id_strings:
                    outputfile.write(line)
                    clusterbonds += 1

            #Angles
            for line in inputlines[2*numatoms+numheaderlines+2*numskips+numbonds:2*numatoms+numheaderlines+3*numskips+numbonds]:
                outputfile.write(line)
            for line in inputlines[2*numatoms+numheaderlines+3*numskips+numbonds:]:
                split = str.split(line)
                if len(split) < 3:
                    pass
                elif split[2] in id_strings and split[3] in id_strings and split[4] in id_strings:
                    outputfile.write(line)
                    clusterangles += 1


            outputfile.close()



            #Correcting header
            if clusteratoms!= len(ids):
                print('This is bad! Length of IDs vector: {} Clusteratoms:{} on '.format(len(ids),clusteratoms) + string)
            sb.call(["gsed","-i",'s:[0-9][0-9][0-9][0-9][0-9]\ atoms:'+"{}".format(clusteratoms)+'\ atoms:',filename+string+'.data'])
            sb.call(["gsed","-i",'s:[0-9][0-9][0-9][0-9][0-9]\ bonds:'+"{}".format(clusterbonds)+'\ bonds:',filename+string+'.data'])
            sb.call(["gsed","-i",'s:[0-9][0-9][0-9][0-9][0-9]\ angles:'+"{}".format(clusterangles)+'\ angles:',filename+string+'.data'])


if __name__ == "__main__":
   usage = '''Usage: interior_exterior filename cutoff_int cutoff_ext
No extension for filename. Cutoffs are optional, must be integer of monos within 2.5sigma and int > ext.
   '''
   
   if len(sys.argv) < 2:
       print(usage)
   elif len(sys.argv) < 3:
       filename = sys.argv[1]
       criteria_int = int(60)
       criteria_ext = int(12)
       calculate_int_ext(filename,criteria_int,criteria_ext)
   elif len(sys.argv) < 4:
       filename = sys.argv[1]
       criteria_int = int(sys.argv[2])
       criteria_ext = int(12)
       calculate_int_ext(filename,criteria_int,criteria_ext)
   else: 
       filename = sys.argv[1]
       criteria_int = int(sys.argv[2])
       criteria_ext = int(sys.argv[3])
       calculate_int_ext(filename,criteria_int,criteria_ext)
