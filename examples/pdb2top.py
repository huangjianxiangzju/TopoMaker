#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 14 23:21:53 2022

@author: huangjianxiang
"""

import MDAnalysis as mda
import sys
import numpy as np

#The CONECT info in the pdb is important!
#It can be read by the mda package
filename="11H_ter_better5.pdb"

u=mda.Universe(filename)

all=u.select_atoms("all")

#all.bonds.
bonds=all.bonds
conect=bonds.dump_contents()+1

num,_=np.shape(conect)

#generate a tcl file to check the conection using VMD
data=open("check.tcl","w")
for i in range(num):
    a=conect[i][0]
    b=conect[i][1]
    data.write("topo addbond "+str(a-1)+"\t"+str(b-1)+"\n")
data.close

#CHECKING CONECTION
for i in range(num):
    a=conect[i][0]
    b=conect[i][1]
    aa=u.select_atoms("bynum "+str(a))
    bb=u.select_atoms("bynum "+str(b))
#    print(aa.names[0],bb.names[0])
    if "O" in aa.names[0] and "O" in bb.names[0]: 
        print(a,b)

#ANGLE SECTION
anglelist=[]
for i in range(num):
    a=conect[i][0]
    b=conect[i][1]
    for j in range(len(conect)):
        c=conect[j][0]
        d=conect[j][1]
        if ( a ==c ) and (b != c and b!=d):
            if [d,a,b] not in anglelist:
                anglelist.append([b,a,d])
        elif (a==d) and (b != c and b!=d):
            if [c,a,b] not in anglelist:
                anglelist.append([b,a,c])
        elif ( b ==c ) and (a != c and a!=d):
            if [d,b,a] not in anglelist:
                anglelist.append([a,b,d])
        elif (b==d) and (a != c and a!=d):
            if [c,b,a] not in anglelist:
                anglelist.append([a,b,c])
#DIHE SECTION
dihelist=[]
for i in range(len(anglelist)):
    a=anglelist[i][0]
    b=anglelist[i][1]
    c=anglelist[i][2]
    
    for j in range(i,len(anglelist)):
        b1=anglelist[j][0]
        c1=anglelist[j][1]
        d1=anglelist[j][2]
        if b==b1 and c==c1:
            dihelist.append([a,b,c,d1])

#--------------------
atom_head="""
[ moleculetype ]
; molname  nrexcl
CSC	    3
[ atoms ]
;  nr    type       resnr  residu    atom    cgnr        charge  
"""

"""
;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
     1   Si     1   TMA   SI1    1     1.362350     28.09000 ; qtot 1.362
     2   do     1   TMA    O1    2     0.000000     16.00000 ; qtot 0.497
     3   oh     1   TMA    O2    3    -0.865839     16.00000 ; qtot -0.369
     4   oh     1   TMA    O3    4    -0.865839     16.00000 ; qtot -1.235
     5   c3     1   TMA    C1    5    -0.430889     12.01000 ; qtot -1.666
     6   c3     1   TMA    C2    6     0.273215     12.01000 ; qtot -1.393
     7   c3     1   TMA    C3    7    -0.197569     12.01000 ; qtot -1.590
     8   n4     1   TMA    N1    8     0.141050     14.01000 ; qtot -1.449
     9   c3     1   TMA    C4    9    -0.369733     12.01000 ; qtot -1.819
    10   c3     1   TMA    C5   10    -0.369733     12.01000 ; qtot -2.189
    11   c3     1   TMA    C6   11    -0.369733     12.01000 ; qtot -2.559
    12   hc     1   TMA    H1   12    -0.004074      1.00800 ; qtot -2.563
    13   hc     1   TMA    H2   13    -0.004074      1.00800 ; qtot -2.567
    14   hx     1   TMA    H3   14     0.131751      1.00800 ; qtot -2.435
    15   hx     1   TMA    H4   15     0.131751      1.00800 ; qtot -2.303
    16   hx     1   TMA    H5   16     0.186220      1.00800 ; qtot -2.117
    17   hx     1   TMA    H6   17     0.186220      1.00800 ; qtot -1.931
    18   hx     1   TMA    H7   18     0.186220      1.00800 ; qtot -1.745
    19   hx     1   TMA    H8   19     0.186220      1.00800 ; qtot -1.558
    20   hx     1   TMA    H9   20     0.186220      1.00800 ; qtot -1.372
    21   hx     1   TMA   H10   21     0.186220      1.00800 ; qtot -1.186
    22   hx     1   TMA   H11   22     0.186220      1.00800 ; qtot -1.000
    23   hx     1   TMA   H12   23     0.186220      1.00800 ; qtot -0.813
    24   hx     1   TMA   H13   24     0.186220      1.00800 ; qtot -0.627
    25   dh     1   TMA   H14   25     0.000000      1.00800 ; qtot -0.141
    26   ho     1   TMA   H15   26     0.486133      1.00800 ; qtot 0.345
    27   ho     1   TMA   H16   27     0.486133      1.00800 ; qtot 0.831
    28   hc     1   TMA   H17   28     0.084412      1.00800 ; qtot 0.916
    29   hc     1   TMA   H18   29     0.084412      1.00800 ; qtot 1.000
"""

"""
   260   Si    12   TEO   SI1  260     1.649764     28.09000 ; qtot 9.650
   261   oh    12   TEO    O1  261    -0.909629     16.00000 ; qtot 8.740
   262   oh    12   TEO    O2  262    -0.909629     16.00000 ; qtot 7.831
   263   oh    12   TEO    O3  263    -0.909629     16.00000 ; qtot 6.921
   264   oh    12   TEO    O4  264    -0.909629     16.00000 ; qtot 6.011
   265   ho    12   TEO    H1  265     0.497187      1.00800 ; qtot 6.508
   266   ho    12   TEO    H2  266     0.497187      1.00800 ; qtot 7.006
   267   ho    12   TEO    H3  267     0.497187      1.00800 ; qtot 7.503
   268   ho    12   TEO    H4  268     0.497187      1.00800 ; qtot 8.000

"""

atom_section=[]
for i in range(all.n_atoms):
    if all.resnames[i]=='TEO':
        if "SI1" in all.names[i]:
            atom_section.append(str(i+1)+"\t"+"Si"+"\t"+str(i+1)+"\t"+"TEO"+"\t"+"SI1"+"\t"+str(i+1)+"\t1.37011\t28.09000\n")
        elif "O" in all.names[i]:
            atom_section.append(str(i+1)+"\t"+"oh"+"\t"+str(i+1)+"\t"+"TEO"+"\t"+"O"+"\t"+str(i+1)+"\t-0.50000\t16.00000\n")
        elif "H" in all.names[i]:
            atom_section.append(str(i+1)+"\t"+"ho"+"\t"+str(i+1)+"\t"+"TEO"+"\t"+"H"+"\t"+str(i+1)+"\t0.18000\t1.00800\n")
        else:
            print("Warning:TEO")
    else:
        if "SI1" in all.names[i]:
            atom_section.append(str(i+1)+"\t"+"Si"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"SI1"+"\t"+str(i+1)+"\t1.37011\t28.09000\n")
        elif "O" in all.names[i]:
            atom_section.append(str(i+1)+"\t"+"oh"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"O"+"\t"+str(i+1)+"\t-0.50000\t16.00000\n")
        elif "C1" in all.names[i]:
            atom_section.append(str(i+1)+"\t"+"c3"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"C1"+"\t"+str(i+1)+"\t-0.27000\t12.01000\n")
        elif "C2" in all.names[i]:
            atom_section.append(str(i+1)+"\t"+"c3"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"C2"+"\t"+str(i+1)+"\t0.25000\t12.01000\n")
        elif "C3" in all.names[i]:
            atom_section.append(str(i+1)+"\t"+"c3"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"C3"+"\t"+str(i+1)+"\t-0.18000\t12.01000\n")
        elif "C4" in all.names[i]:
            atom_section.append(str(i+1)+"\t"+"c3"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"C4"+"\t"+str(i+1)+"\t-0.35000\t12.01000\n")
        elif "C5" in all.names[i]:
            atom_section.append(str(i+1)+"\t"+"c3"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"C5"+"\t"+str(i+1)+"\t-0.35000\t12.01000\n")
        elif "C6" in all.names[i]:
            atom_section.append(str(i+1)+"\t"+"c3"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"C6"+"\t"+str(i+1)+"\t-0.35000\t12.01000\n")
        elif "H1" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hc"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H1"+"\t"+str(i+1)+"\t-0.0041\t1.00800\n")
        elif "H2" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hc"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H2"+"\t"+str(i+1)+"\t-0.0041\t1.00800\n")
        elif "H3" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hx"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H3"+"\t"+str(i+1)+"\t0.13000\t1.00800\n")
        elif "H4" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hx"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H4"+"\t"+str(i+1)+"\t0.13000\t1.00800\n")
        elif "H5" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hx"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H5"+"\t"+str(i+1)+"\t0.18620\t1.00800\n")
        elif "H6" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hx"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H6"+"\t"+str(i+1)+"\t0.18620\t1.00800\n")
        elif "H7" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hx"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H7"+"\t"+str(i+1)+"\t0.186220\t1.00800\n")
        elif "H8" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hx"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H8"+"\t"+str(i+1)+"\t0.186220\t1.00800\n")
        elif "H9" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hx"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H9"+"\t"+str(i+1)+"\t0.186220\t1.00800\n")
        elif "H10" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hx"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H10"+"\t"+str(i+1)+"\t0.186220\t1.00800\n")
        elif "H11" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hx"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H11"+"\t"+str(i+1)+"\t0.186220\t1.00800\n")
        elif "H12" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hx"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H12"+"\t"+str(i+1)+"\t0.186220\t1.00800\n")
        elif "H13" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hx"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H13"+"\t"+str(i+1)+"\t0.186220\t1.00800\n")
        elif "H14" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"ho"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H14"+"\t"+str(i+1)+"\t0.18000\t1.00800\n")
        elif "H15" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"ho"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H15"+"\t"+str(i+1)+"\t0.18000\t1.00800\n")
        elif "H16" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"ho"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H16"+"\t"+str(i+1)+"\t0.18000\t1.00800\n")
        elif "H17" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hc"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H17"+"\t"+str(i+1)+"\t0.08400\t1.00800\n")
        elif "H18" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"hc"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"H18"+"\t"+str(i+1)+"\t0.08400\t1.00800\n")
        elif "N1" == all.names[i]:
            atom_section.append(str(i+1)+"\t"+"n4"+"\t"+str(i+1)+"\t"+"TMA"+"\t"+"N1"+"\t"+str(i+1)+"\t0.143000\t14.01000\n")
        else:
            print("Warning:TMA")
            
 
bond_head="""
 
[ bonds ]
;  ai  aj   funct
"""

pair_head="""

[ pairs ]
;   ai     aj    funct
"""


angle_head="""
[ angles ]
;   ai     aj     ak    funct   theta         cth
"""

dihe_head="""
[ dihedrals ] ; propers
; for gromacs 4.5 or higher, using funct 9
;    i      j      k      l   func   phase     kd      pn
"""

bond_section=[]
for i in range(len(conect)):
    a=conect[i][0]
    b=conect[i][1]
    bond_section.append(str(conect[i][0])+"\t"+str(conect[i][1])+"\t1\n")

angle_section=[]
for i in range(len(anglelist)):
    a=anglelist[i][0]
    b=anglelist[i][1]
    c=anglelist[i][2]
    
    angle_section.append(str(anglelist[i][0])+"\t"+str(anglelist[i][1])+"\t"+str(anglelist[i][2])+"\t1\n")
        
dihe_section=[]
pair_section=[]
for i in range(len(dihelist)):
    a=dihelist[i][0]
    b=dihelist[i][1]
    c=dihelist[i][2]
    d=dihelist[i][3]
    
    dihe_section.append(str(dihelist[i][0])+"\t"+str(dihelist[i][1])+"\t"+str(dihelist[i][2])+"\t"+str(dihelist[i][3])+"\t9\n")
    pair_section.append(str(dihelist[i][0])+"\t"+str(dihelist[i][3])+"\t1\n")

        
output=open("output.itp","w")
output.writelines(atom_head)
output.writelines(atom_section)
output.writelines(bond_head)
output.writelines(bond_section)
output.writelines(pair_head)
output.writelines(pair_section)
output.writelines(angle_head)
output.writelines(angle_section)
output.writelines(dihe_head)
output.writelines(dihe_section)
output.close()

##################MASS-CHARGE-REFINE--->TO DO!
