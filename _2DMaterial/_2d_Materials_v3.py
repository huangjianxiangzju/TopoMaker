# -*- coding: utf-8 -*-
"""
Created on Sat Jan  5 17:57:04 2019

@author: huangjianxiang

huangjianxiang@zju.edu.cn
MoS2 , graphene , C2N structures can be generated.

I add functions to generate Ti2C and Ti3C2.

I add functions to generate h-BN


beta version. 
"""
import math
import numpy as np
import matplotlib.pyplot as plt

def mos2_diamond(m):
    a=3.150# 边长，单位angstrom
    n=m+1
    output=open("mos2_output_diamond"+str(m)+".xyz",'w')
    z1=1.488301;z2=3.075000;z3=4.661700#z1 z3 内层和外层S的高度；z2 中间Mo的高度
    coord1=np.zeros((n**2,2)) #the sulfor atoms
    coord2=np.zeros((n**2,2)) #the Mo atoms
    for i in range(0,n):
        for j in range(0,n):
            if i==0:
                if j==0:
                    coord1[i*n+j,0]=0
                    coord1[i*n+j,1]=0
                    coord2[i*n+j,0]=0.5*a
                    coord2[i*n+j,1]=0.5*a/np.sqrt(3)
                else:
                    coord1[i*n+j,0]=coord1[i*n+j-1,0]+a
                    coord1[i*n+j,1]=coord1[i*n+j-1,1]
                    coord2[i*n+j,0]=coord2[i*n+j-1,0]+a
                    coord2[i*n+j,1]=coord2[i*n+j-1,1]
            else:
                if j==0:
                    coord1[i*n+j,0]=coord1[(i-1)*n+j,0]+0.5*a
                    coord1[i*n+j,1]=coord1[(i-1)*n+j,1]+np.sqrt(3)/2*a
                    coord2[i*n+j,0]=coord2[(i-1)*n+j,0]+a*0.5
                    coord2[i*n+j,1]=coord2[(i-1)*n+j,1]+np.sqrt(3)*a*0.5
                else:
                    coord1[i*n+j,0]=coord1[i*n+j-1,0]+a
                    coord1[i*n+j,1]=coord1[i*n+j-1,1]
                    coord2[i*n+j,0]=coord2[i*n+j-1,0]+a
                    coord2[i*n+j,1]=coord2[i*n+j-1,1]
    plt.scatter(coord1[:,0],coord1[:,1])
    plt.scatter(coord2[:,0],coord2[:,1])
    output.write(str(3*n**2))
    output.write("\n\n")
    for i in range(n**2):
        output.write( " S       "+'%6.5f'%coord1[i,0] +"\t "+ '%6.5f'%coord1[i,1] +"\t "+str(z1)+"\n")
        output.write( "Mo       "+'%6.5f'%coord2[i,0] +"\t "+ '%6.5f'%coord2[i,1] +"\t "+str(z2)+"\n")
        output.write( " S       "+'%6.5f'%coord1[i,0] +"\t "+ '%6.5f'%coord1[i,1] +"\t "+str(z3)+"\n")
    output.close()    
 
def mos2(n,m):
    num=m*n*12  #the total number of atoms
    a=6.3  # x axis in angstrom
    b=5.45596  # y axis in angstrom
    output1=open("mos2_output_"+str(n)+"_"+str(m)+".gro",'w')
    coord=np.array([[15.750000,16.367880, 1.488301],  
             [17.325001,17.277210, 3.075000],
             [15.750000,16.367880, 4.661700],
             [18.900000,16.367880, 1.488301],
             [20.475000,17.277210, 3.075000],
             [18.900000,16.367880, 4.661700],
             [15.750000,20.005190, 3.075000],
             [17.325001,19.095860, 1.488301],
             [18.900000,20.005190, 3.075000],
             [17.325001,19.095860, 4.661700],
             [20.475000,19.095860, 1.488301],
             [20.475000,19.095860, 4.661700]]) #square unit 
    unit_element=[" S","Mo"," S"," S","Mo"," S","Mo",\
                  " S","Mo"," S"," S"," S"] #the atom name of the square unit
    coord[:,0]=coord[:,0]  -  15*np.ones_like(coord[:,0])
    coord[:,1]=coord[:,1]  -  16*np.ones_like(coord[:,1])
    coord[:,2]=coord[:,2]  +  2*np.ones_like(coord[:,2])  
    #change the coordinate to fit the PBC
    coord1=np.zeros((12,3,m))
    coord2=np.zeros((12,3,m,n))
    for i in range(m):  #x direction
        coord1[:,0,i]=coord[:,0]+i*a*np.ones_like(coord1[:,0,i])
        coord1[:,1,i]=coord[:,1]
        coord1[:,2,i]=coord[:,2]
    for j in range(n):  #y direction
        coord2[:,1,:,j]=coord1[:,1,:]+j*b*np.ones_like(coord2[:,1,:,j])
        coord2[:,0,:,j]=coord1[:,0,:]
        coord2[:,2,:,j]=coord1[:,2,:]
    output1.write("Mo-S bond length=3.150  huangjianxiang\n")
    output1.write(str(num)+"\n")
    #below PBC
    y=(n)*5.45596 #in angstrom
    x=(m)*6.3    #in angstrom
    z=10   #in angstrom
    print("consider using the PBC xyz "+'%6.5f'%y+"   "+'%6.5f'%x+"   "+'%6.5f'%z+" in angstrom")
    id=1
    for i in range(m):
        for j in range(n):
            for k in range(12):
                output1.write( "    mos2   "+unit_element[k]+"_2"+'%5d'%id+'%8.4f'%(coord2[k,0,i,j]*0.1) + '%8.4f'%(coord2[k,1,i,j]*0.1) +'%8.4f'%(coord2[k,2,i,j]*0.1)+"\n")
                id+=1
                if id ==100000:
                    id =1 
    output1.write("  "+'%10.6f'%(x*0.1)+""+'%12.6f'%(y*0.1)+""+'%12.6f'%(z*0.1)+"\n")
    output1.close()      

def gra(m,n):
    num=m*n*4 #build a square shape first
    a=1.418 #angstrom
    output1=open("gra_output_"+str(m)+"_"+str(n)+".gro",'w')
    coord=np.zeros((4,2))
    coord1=np.zeros((n*4,2))
    coord2=np.zeros((n*4,2,m))
    coord[0,:]=0.00,0.00
    coord[1,:]=1.228024,0.709
    coord[2,:]=1.228024,2.127
    coord[3,:]=0.00,2.836
    for i in range(n):
        coord1[4*i:4*(i+1),0]=coord[:,0]+i*a*np.sqrt(3)*np.ones_like(coord[:,0])
        coord1[4*i:4*(i+1),1]=coord[:,1]
    for j in range(m):
        coord2[:,0,j]=coord1[:,0]
        coord2[:,1,j]=coord1[:,1]+3*a*j*np.ones_like(coord1[:,1])
    #plt.scatter(coord2[:,0],coord2[:,1])
    output1.write("bond length=1.4  huangjianxiang\n")
    output1.write(str(num)+"\n")
    x=n*0.2456047*10    
    y=m*0.425400 *10 
    z=0.280000*10
    print("consider using the PBC xyz "+'%6.5f'%x+"   "+'%6.5f'%y+"   "+'%6.5f'%z+" in angstrom")
    #output.write("consider using the PBC xyz "+'%6.5f'%x+"   "+'%6.5f'%y+"   "+'%6.5f'%z+"\n")
    id=1
    for i in range(m):
        for j in range(n):
            output1.write( "     gra   "+"Cgra"+'%5d'%id+'%8.4f'%(coord2[j*4-4,0,i]*0.1) + '%8.4f'%(coord2[j*4-4,1,i]*0.1) +'%8.4f'%(1.4000*0.1)+"\n")
            output1.write( "     gra   "+"Cgra"+'%5d'%id+'%8.4f'%(coord2[j*4-3,0,i]*0.1) + '%8.4f'%(coord2[j*4-3,1,i]*0.1) +'%8.4f'%(1.4000*0.1)+"\n")
            output1.write( "     gra   "+"Cgra"+'%5d'%id+'%8.4f'%(coord2[j*4-2,0,i]*0.1) + '%8.4f'%(coord2[j*4-2,1,i]*0.1) +'%8.4f'%(1.4000*0.1)+"\n")
            output1.write( "     gra   "+"Cgra"+'%5d'%id+'%8.4f'%(coord2[j*4-1,0,i]*0.1) + '%8.4f'%(coord2[j*4-1,1,i]*0.1) +'%8.4f'%(1.4000*0.1)+"\n")
            id+=1
            if id ==100000:
                id =1	
    output1.write("  "+'%10.6f'%(x*0.1)+""+'%12.6f'%(y*0.1)+""+'%12.6f'%(z*0.1)+"\n")
    output1.close() 

    output2=open("gra_output_"+str(m)+"_"+str(n)+".itp",'w')
    #HEAD 
    output2.write("[ moleculetype ]\n")
    output2.write("; name  nrexcl\n")
    output2.write("gra 3\n")
    #ATOM
    output2.write("\n[ atoms ]\n;   nr    type   resnr  residu    atom    cgnr        charge          mass\n")  
    unit_element=["C","C","C","C"]
    charge=[0.0,0.0,0.0,0.0]
    id=0
    for i in range(m):
        for j in range(n):
            for k in range(4):
                id+=1
                output2.write( '%5d'%id+"   "+unit_element[k]+"  1    gra   " +unit_element[k].lower() +"\t1\t"+str(charge[k])+"\n")
                if id ==100000:
                    id =1
    output2.close()

#function: add h-BN 
def hbn(m,n):
    num=m*n*4 #build a square shape first
    a=1.451 #angstrom
    output1=open("hbn_output_"+str(m)+"_"+str(n)+".gro",'w')
    coord=np.zeros((4,2))
    coord1=np.zeros((n*4,2))
    coord2=np.zeros((n*4,2,m))
    coord[0,:]=0.00,0.00  #B
    coord[1,:]=  math.sqrt(3)*0.5*a ,0.5*a #N
    coord[2,:]=  math.sqrt(3)*0.5*a ,1.5*a #B
    coord[3,:]=0.00,2*a  #N
    for i in range(n):
        coord1[4*i:4*(i+1),0]=coord[:,0]+i*a*np.sqrt(3)*np.ones_like(coord[:,0])
        coord1[4*i:4*(i+1),1]=coord[:,1]
    for j in range(m):
        coord2[:,0,j]=coord1[:,0]
        coord2[:,1,j]=coord1[:,1]+3*a*j*np.ones_like(coord1[:,1])
    #plt.scatter(coord2[:,0],coord2[:,1])
    output1.write("bond length=1.4  huangjianxiang\n")
    output1.write(str(num)+"\n")
    x=n*0.2456047*10    
    y=m*0.425400 *10 
    z=0.280000*10
    print("consider using the PBC xyz "+'%6.5f'%x+"   "+'%6.5f'%y+"   "+'%6.5f'%z+" in angstrom")
    #output.write("consider using the PBC xyz "+'%6.5f'%x+"   "+'%6.5f'%y+"   "+'%6.5f'%z+"\n")
    id=1
    for i in range(m):
        for j in range(n):
            output1.write( "     hbn   "+"Bhbn"+'%5d'%id+'%8.4f'%(coord2[j*4-4,0,i]*0.1) + '%8.4f'%(coord2[j*4-4,1,i]*0.1) +'%8.4f'%(1.4000*0.1)+"\n")
            output1.write( "     hbn   "+"Nhbn"+'%5d'%id+'%8.4f'%(coord2[j*4-3,0,i]*0.1) + '%8.4f'%(coord2[j*4-3,1,i]*0.1) +'%8.4f'%(1.4000*0.1)+"\n")
            output1.write( "     hbn   "+"Bhbn"+'%5d'%id+'%8.4f'%(coord2[j*4-2,0,i]*0.1) + '%8.4f'%(coord2[j*4-2,1,i]*0.1) +'%8.4f'%(1.4000*0.1)+"\n")
            output1.write( "     hbn   "+"Nhbn"+'%5d'%id+'%8.4f'%(coord2[j*4-1,0,i]*0.1) + '%8.4f'%(coord2[j*4-1,1,i]*0.1) +'%8.4f'%(1.4000*0.1)+"\n")
            id+=1
            if id ==100000:
                id =1   
    output1.write("  "+'%10.6f'%(x*0.1)+""+'%12.6f'%(y*0.1)+""+'%12.6f'%(z*0.1)+"\n")
    output1.close()   

    output2=open("hbn_output_"+str(m)+"_"+str(n)+".itp",'w')
    #HEAD 
    output2.write("[ moleculetype ]\n")
    output2.write("; name  nrexcl\n")
    output2.write("hbn 3\n")
    #ATOM
    output2.write("\n[ atoms ]\n;   nr    type   resnr  residu    atom    cgnr        charge          mass\n")
    unit_element=["B","N","B","N"]
    charge=[0.4,-0.4,0.4,-0.4]
    id=0
    for i in range(m):
        for j in range(n):
            for k in range(4):
                id+=1
                output2.write( '%5d'%id+"   "+unit_element[k]+"  1    hbn   " +unit_element[k].lower() +"\t1\t"+str(charge[k])+"\n")
                if id ==100000:
                    id =1
    output2.close()
                #if coord2[k,2,i,j]==6.159:
                #    id+=1
                #    output2.write( '%5d'%id+"    O_2  1    hbn    o_2\t1\t-0.79\n")
                #    id+=1
                #    output2.write( '%5d'%id+"    H_2  1    hbn    h_2\t1\t 0.35\n")
                #elif coord2[k,2,i,j]==1.409:
                #    id+=1
                #    output2.write( '%5d'%id+"    O_2  1    hbn    o_2\t1\t-0.79\n")
                #    id+=1
                #    output2.write( '%5d'%id+"    H_2  1    hbn    h_2\t1\t 0.35\n")
#
def c2n(m,n):
    num=m*n*36
    #output=open("c2n_output_"+str(m)+"_"+str(n)+".xyz",'w')
    output1=open("c2n_output_"+str(m)+"_"+str(n)+".gro",'w')
    a=6*1.4  # wideth
    b=6*np.sqrt(3)*1.4  # length
    unit=np.zeros((36,2))
    coord=np.zeros((36,2,m,n))
    unit[0,:]   =    13.94436 , 2.80000   
    unit[1,:]   =    12.73436 , 3.50000   
    unit[2,:]   =     1.82000 , 2.80000   
    unit[3,:]   =     0.61000 , 3.50000   
    unit[4,:]   =     4.24487 , 2.80000   
    unit[5,:]   =     3.03487 , 3.50000   
    unit[6,:]   =     5.45974 , 3.50000   
    unit[7,:]   =    11.51948 , 2.80000   
    unit[8,:]   =    10.30948 , 3.50000   
    unit[9,:]   =    12.73436 , 4.90000   
    unit[10,:]  =    12.73436 , 7.70000   
    unit[11,:]  =     3.03487 , 4.90000   
    unit[12,:]  =     4.24487 , 5.60000   
    unit[13,:]  =     4.24487 , 7.00000   
    unit[14,:]  =     3.03487 , 7.70000   
    unit[15,:]  =     5.45974 , 4.90000   
    unit[16,:]  =     6.66974 , 5.60000   
    unit[17,:]  =     6.66974 , 7.00000   
    unit[18,:]  =     5.45974 , 7.70000   
    unit[19,:]  =     7.88461 , 4.90000   
    unit[20,:]  =     9.09461 , 5.60000   
    unit[21,:]  =     9.09461 , 7.00000   
    unit[22,:]  =     7.88461 , 7.70000   
    unit[23,:]  =    10.30948 , 4.90000   
    unit[24,:]  =    11.51948 , 5.60000   
    unit[25,:]  =    11.51948 , 7.00000   
    unit[26,:]  =    10.30948 , 7.70000   
    unit[27,:]  =    12.73436 , 9.10000   
    unit[28,:]  =    13.94436 , 9.80000   
    unit[29,:]  =     0.61000 , 9.10000   
    unit[30,:]  =     1.82000 , 9.80000   
    unit[31,:]  =     3.03487 , 9.10000   
    unit[32,:]  =     4.24487 , 9.80000   
    unit[33,:]  =     5.45974 , 9.10000   
    unit[34,:]  =    10.30948 , 9.10000   
    unit[35,:]  =    11.51948 , 9.80000  
    unit_element=["C","C","C","N","C","C","N","C","N","N","N","N","C","C","N",\
                  "C","C","C","C","N","C","C","N","C","C","C","C","C","C","N",\
                  "C","C","C","N","N","C"]
    right=np.zeros_like(unit)
    up=np.zeros_like(unit)
    right[:,0]=b*np.ones_like(right[:,0])
    up[:,1]=a*np.ones_like(right[:,1])
    for i in range(n):
        for j in range(m):
            if i==0 and j==0:
                coord[:,:,j,i]=unit
            else:
                coord[:,:,j,i]=unit+i*right+j*up
    #output.write(str(num)+"\n")  输出xyz格式已经过时，不再需要
    #output.write("\n")
	#above xyz
	#below gro
    output1.write("bond length=1.4  huangjianxiang\n")
    output1.write(str(num)+"\n")
    x=n*1.454922*10    
    y=m*0.840000 *10 
    z=0.280000*10
    print("consider using the PBC xyz "+'%6.5f'%x+"   "+'%6.5f'%y+"   "+'%6.5f'%z+" in angstrom")
    id=1
    for i in range(n):
        for j in range(m):
            for k in range(36):
                #output.write( unit_element[k]+ "/t"+'%6.5f'%coord[k,0,j,i] +"\t "+ '%6.5f'%(coord[k,1,j,i]-2.1) +"\t 0.0000\n")
                output1.write( "     c2n   "+unit_element[k]+"c2n"+'%5d'%id+'%8.4f'%(coord[k,0,j,i]*0.1) + '%8.4f'%((coord[k,1,j,i]-2.1)*0.1) +'%8.4f'%(1.4000*0.1)+"\n")
                id+=1
                if id ==100000:
                    id =1				
    #output.close()
    output1.write("  "+'%10.6f'%(x*0.1)+""+'%12.6f'%(y*0.1)+""+'%12.6f'%(z*0.1)+"\n")
    output1.close()  	

def ti2c(m,n):
    num=m*n*6
    output1=open("ti2c_output_"+str(m)+"_"+str(n)+".gro",'w')
    a=3.08  # wideth unit:angstrom
    b=5.34  # length
    unit=np.zeros((6,3))
    coord=np.zeros((6,3,m,n))
    coord=np.array([[0.835,   0.436,   1.260 ],  
             [0.835  , 2.216  , 3.550],
             [2.375  , 1.326  , 2.410],
             [2.375  , 3.116  , 1.260],
             [0.835  , 4.006  , 2.410],
             [2.375  , 4.896  , 3.550]]) #square unit
    unit_element=["Ti","Ti"," C","Ti"," C","Ti"]
    coord1=np.zeros((6,3,m))
    coord2=np.zeros((6,3,m,n))
    for i in range(m):  #x direction
        coord1[:,0,i]=coord[:,0]+i*a*np.ones_like(coord1[:,0,i])
        coord1[:,1,i]=coord[:,1]
        coord1[:,2,i]=coord[:,2]
    for j in range(n):  #y direction
        coord2[:,1,:,j]=coord1[:,1,:]+j*b*np.ones_like(coord2[:,1,:,j])
        coord2[:,0,:,j]=coord1[:,0,:]
        coord2[:,2,:,j]=coord1[:,2,:]
    output1.write("Ti-C bond length=2.117  huangjianxiang\n")
    output1.write(str(num)+"\n")
    #below PBC
    y=(n)*b #in angstrom
    x=(m)*a    #in angstrom
    z= 4.814   #in angstrom
    print("consider using the PBC xyz "+'%6.5f'%y+"   "+'%6.5f'%x+"   "+'%6.5f'%z+" in angstrom")
    id=1
    for i in range(m):
        for j in range(n):
            for k in range(6):
                output1.write( "    ti2c   "+unit_element[k]+"_2"+'%5d'%id+'%8.4f'%(coord2[k,0,i,j]*0.1) + '%8.4f'%(coord2[k,1,i,j]*0.1) +'%8.4f'%(coord2[k,2,i,j]*0.1)+"\n")
                id+=1
                if id ==100000:
                    id =1 
    output1.write("  "+'%10.6f'%(x*0.1)+""+'%12.6f'%(y*0.1)+""+'%12.6f'%(z*0.1)+"\n")
    output1.close()      


def ti3c2(m,n):
    num=m*n*10
    output1=open("ti3c2_output_"+str(m)+"_"+str(n)+".gro",'w')
    a=3.07  # wideth unit:angstrom
    b=5.32  # length
    unit=np.zeros((10,3))
    coord=np.zeros((10,3,m,n))
    coord=np.array([[0.765 ,  0.268,   3.779 ],  
             [2.305  , 1.158  , 1.409],
             [0.765  , 2.038  , 6.159],
             [0.765  , 2.038  , 2.489],
             [2.305  , 1.158  , 5.079],
             [0.765  , 3.808  , 1.409],
             [2.305  , 2.928  , 3.779],
             [0.765  , 3.808  , 5.079],
             [2.305  , 4.698  , 6.159],
             [2.305  , 4.698  , 2.489]]) #square unit
    unit_element=["Ti","Ti","Ti"," C"," C","Ti","Ti"," C","Ti"," C"]
    coord1=np.zeros((10,3,m))
    coord2=np.zeros((10,3,m,n))
    for i in range(m):  #x direction
        coord1[:,0,i]=coord[:,0]+i*a*np.ones_like(coord1[:,0,i])
        coord1[:,1,i]=coord[:,1]
        coord1[:,2,i]=coord[:,2]
    for j in range(n):  #y direction
        coord2[:,1,:,j]=coord1[:,1,:]+j*b*np.ones_like(coord2[:,1,:,j])
        coord2[:,0,:,j]=coord1[:,0,:]
        coord2[:,2,:,j]=coord1[:,2,:]
    output1.write("Ti-C bond length=2.117  huangjianxiang\n")
    output1.write(str(num)+"\n")
    #below PBC
    y=(n)*b #in angstrom
    x=(m)*a    #in angstrom
    z= 7.566  #in angstrom
    print("consider using the PBC xyz "+'%6.5f'%y+"   "+'%6.5f'%x+"   "+'%6.5f'%z+" in angstrom")
    id=1
    for i in range(m):
        for j in range(n):
            for k in range(10):
                output1.write( "    ti2c   "+unit_element[k]+"_2"+'%5d'%id+'%8.4f'%(coord2[k,0,i,j]*0.1) + '%8.4f'%(coord2[k,1,i,j]*0.1) +'%8.4f'%(coord2[k,2,i,j]*0.1)+"\n")
                id+=1
                if id ==100000:
                    id =1 
    output1.write("  "+'%10.6f'%(x*0.1)+""+'%12.6f'%(y*0.1)+""+'%12.6f'%(z*0.1)+"\n")
    output1.close()    


def ti3c2oh2(m,n):
    num=m*n*18
    output1=open("ti3c2oh2_output_"+str(m)+"_"+str(n)+".gro",'w')
    a=3.07  # wideth unit:angstrom
    b=5.32  # length
    unit=np.zeros((10,3))
    coord=np.zeros((10,3,m,n))
    coord=np.array([[0.765 ,  0.268,   3.779 ],  
             [2.305  , 1.158  , 1.409],
             [0.765  , 2.038  , 6.159],
             [0.765  , 2.038  , 2.489],
             [2.305  , 1.158  , 5.079],
             [0.765  , 3.808  , 1.409],
             [2.305  , 2.928  , 3.779],
             [0.765  , 3.808  , 5.079],
             [2.305  , 4.698  , 6.159],
             [2.305  , 4.698  , 2.489]]) #square unit
    unit_element=["ti_i","ti_o","ti_o"," c_2"," c_2","ti_o","ti_i"," c_2","ti_o"," c_2"]
    coord1=np.zeros((10,3,m))
    coord2=np.zeros((10,3,m,n))
    for i in range(m):  #x direction
        coord1[:,0,i]=coord[:,0]+i*a*np.ones_like(coord1[:,0,i])
        coord1[:,1,i]=coord[:,1]
        coord1[:,2,i]=coord[:,2]
    for j in range(n):  #y direction
        coord2[:,1,:,j]=coord1[:,1,:]+j*b*np.ones_like(coord2[:,1,:,j])
        coord2[:,0,:,j]=coord1[:,0,:]
        coord2[:,2,:,j]=coord1[:,2,:]
    output1.write("Ti-C bond length=2.117  huangjianxiang\n")
    output1.write(str(num)+"\n")
    #below PBC
    y=(n)*b #in angstrom
    x=(m)*a    #in angstrom
    z= 7.566  #in angstrom
    print("consider using the PBC xyz "+'%6.5f'%y+"   "+'%6.5f'%x+"   "+'%6.5f'%z+" in angstrom")
    id=0
    for i in range(m):
        for j in range(n):
            for k in range(10):
                id+=1
                output1.write( "    ti2c   "+unit_element[k]+'%5d'%id+'%8.4f'%(coord2[k,0,i,j]*0.1) + '%8.4f'%(coord2[k,1,i,j]*0.1) +'%8.4f'%(coord2[k,2,i,j]*0.1)+"\n")
                if id ==100000:
                    id =1
                if coord2[k,2,i,j]==6.159:
                    id+=1
                    output1.write( "    ti2c   "+" o"+"_2"+'%5d'%id+'%8.4f'%(coord2[k,0,i,j]*0.1) + '%8.4f'%(coord2[k,1,i,j]*0.1) +'%8.4f'%(coord2[k,2,i,j]*0.1+0.2)+"\n")
                    id+=1
                    output1.write( "    ti2c   "+" h"+"_2"+'%5d'%id+'%8.4f'%(coord2[k,0,i,j]*0.1) + '%8.4f'%(coord2[k,1,i,j]*0.1+0.05) +'%8.4f'%(coord2[k,2,i,j]*0.1+0.25)+"\n")
                elif coord2[k,2,i,j]==1.409:
                    id+=1
                    output1.write( "    ti2c   "+" o"+"_2"+'%5d'%id+'%8.4f'%(coord2[k,0,i,j]*0.1) + '%8.4f'%(coord2[k,1,i,j]*0.1) +'%8.4f'%(coord2[k,2,i,j]*0.1-0.2)+"\n")
                    id+=1
                    output1.write( "    ti2c   "+" h"+"_2"+'%5d'%id+'%8.4f'%(coord2[k,0,i,j]*0.1) + '%8.4f'%(coord2[k,1,i,j]*0.1-0.05) +'%8.4f'%(coord2[k,2,i,j]*0.1-0.25)+"\n")
    output1.write("  "+'%10.6f'%(x*0.1)+""+'%12.6f'%(y*0.1)+""+'%12.6f'%(z*0.1)+"\n")
    output1.close()
    #generate the correspongind itp file
    output2=open("ti3c2oh2_output_"+str(m)+"_"+str(n)+".itp",'w')
    #HEAD 
    output2.write("[ moleculetype ]\n")
    output2.write("; name  nrexcl\n")
    output2.write("Ti3C 3\n")
    #ATOM
    output2.write("\n[ atoms ]\n;   nr    type   resnr  residu    atom    cgnr        charge          mass\n")
    unit_element=["Ti_i","Ti_o","Ti_o"," C_2"," C_2","Ti_o","Ti_i"," C_2","Ti_o"," C_2"]
    charge=[0.64,0.88,0.88,-0.76,-0.76,0.88,0.64,-0.76,0.88,-0.76]
    id=0
    for i in range(m):
        for j in range(n):
            for k in range(10):
                id+=1
                output2.write( '%5d'%id+"   "+unit_element[k]+"  1    ti3   " +unit_element[k].lower() +"\t1\t"+str(charge[k])+"\n")
                if id ==100000:
                    id =1
                if coord2[k,2,i,j]==6.159:
                    id+=1
                    output2.write( '%5d'%id+"    O_2  1    ti3    o_2\t1\t-0.79\n")
                    id+=1
                    output2.write( '%5d'%id+"    H_2  1    ti3    h_2\t1\t 0.35\n")
                elif coord2[k,2,i,j]==1.409:
                    id+=1
                    output2.write( '%5d'%id+"    O_2  1    ti3    o_2\t1\t-0.79\n")
                    id+=1
                    output2.write( '%5d'%id+"    H_2  1    ti3    h_2\t1\t 0.35\n")
    output2.write("[ bonds ]\n;  ai    aj funct           c0           c1\n")
    id=0
    for i in range(m):
        for j in range(n):
            for k in range(10):
                id+=1
                if id ==100000:
                    id =1
                if coord2[k,2,i,j]==6.159:
                    output2.write(str(id)+"\t"+str(id+1)+"\t"+"1\n")
                    id+=1
                    output2.write(str(id)+"\t"+str(id+1)+"\t"+"1\n")
                    id+=1
                elif coord2[k,2,i,j]==1.409:
                    output2.write(str(id)+"\t"+str(id+1)+"\t"+"1\n")
                    id+=1
                    output2.write(str(id)+"\t"+str(id+1)+"\t"+"1\n")
                    id+=1
    output2.write("[ angles ]\n;  ai    aj    ak funct           c0        c1\n")
    id=0
    for i in range(m):
        for j in range(n):
            for k in range(10):
                id+=1
                if id ==100000:
                    id =1
                if coord2[k,2,i,j]==6.159:
                    id+=1
                    output2.write(str(id-1)+"\t"+str(id)+"\t"+str(id+1)+"\t"+"1\n")
                    id+=1
                elif coord2[k,2,i,j]==1.409:
                    id+=1
                    output2.write(str(id-1)+"\t"+str(id)+"\t"+str(id+1)+"\t"+"1\n")
                    id+=1
    output2.close()


