#!/usr/bin/env python
# coding: utf-8

conect=[]
for i in range(24):
    conect.append([1+i*11,2+i*11])
    conect.append([2+i*11,3+i*11])
    conect.append([3+i*11,5+i*11])
    conect.append([2+i*11,4+i*11])
    conect.append([5+i*11,6+i*11])
    conect.append([6+i*11,7+i*11])
    conect.append([6+i*11,8+i*11])
    conect.append([8+i*11,9+i*11])
    conect.append([8+i*11,10+i*11])
    conect.append([4+i*11,11+i*11])
    if i<23:
        conect.append([7+i*11,12+i*11])

anglelist=[]
for i in range(262):
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

"""
PEI (sequence: C1 C2 C3 C4 C5 C6 C7 C8 C9 CX C0)
      C0      C9 CX
      |        \/
      C4       C8
      |        |  
-[-C1-C2-C3-C5-C6-C7-]-

protonated sites
3,7,9,10,11
"""

#left type
ltype=["pq","t","sq","s","s","t","sq","t","pq","pq","pq"]
#main type
mtype=["s","t","sq","s","s","t","sq","t","pq","pq","pq"]
#right type
rtype=["s","t","sq","s","s","t","pq","t","pq","pq","pq"]

#all type
#11*24=264 beads number
pei_atype=[]
for i in range(24):
    if i==0:
        for j in range(len(ltype)):
            pei_atype.append(ltype[j])
    elif i==23:
        for k in range(len(rtype)):
            pei_atype.append(rtype[k])
    else:
        for l in range(len(mtype)):
            pei_atype.append(mtype[l])

"""
24% PEI (sequence: C1 C2 C3 C4 C5 C6 C7 C8 C9 CX C0)
      C0      C9 CX
      |        \/
      C4       C8
      |        |  
-[-C1-C2-C3-C5-C6-C7-]-

protonated sites
3,7,9,11
"""

#left type
ltype=["p","t","sq","s","s","t","sq","t","pq","p","pq"]
#main type
mtype=["s","t","sq","s","s","t","sq","t","pq","p","pq"]
#right type
rtype=["s","t","sq","s","s","t","pq","t","pq","p","pq"]

#all type
#11*24=264 beads number
pei24_atype=[]
for i in range(24):
    if i==0:
        for j in range(len(ltype)):
            pei24_atype.append(ltype[j])
    elif i==23:
        for k in range(len(rtype)):
            pei24_atype.append(rtype[k])
    else:
        for l in range(len(mtype)):
            pei24_atype.append(mtype[l])


"""
87% PEI (sequence: C1 C2 C3 C4 C5 C6 C7 C8 C9 CX C0)
      C0      C9 CX
      |        \/
      C4       C8
      |        |  
-[-C1-C2-C3-C5-C6-C7-]-

protonated sites
3
"""

#left type
ltype=["p","t","sq","s","s","t","s","t","p","p","p"]
#main type
mtype=["s","t","sq","s","s","t","s","t","p","p","p"]
#right type
rtype=["s","t","sq","s","s","t","p","t","p","p","p"]

#all type
#11*24=264 beads number
pei87_atype=[]
for i in range(24):
    if i==0:
        for j in range(len(ltype)):
            pei87_atype.append(ltype[j])
    elif i==23:
        for k in range(len(rtype)):
            pei87_atype.append(rtype[k])
    else:
        for l in range(len(mtype)):
            pei87_atype.append(mtype[l])





info=open("Bond_length.txt")

p2_list=["t","s","p"]
qd_list=["tq","sq","pq"]

atom_head="""
[ moleculetype ]
; molname  nrexcl
PEI	    3

[ atoms ]
;  nr    type       resnr  residu    atom    cgnr        charge  
"""

atom_section=[]
for i in range(24):
    for j in range(11):
        if pei_atype[i*11+j] in p2_list:
            atom_section.append(str(i*11+j+1)+"\t"+"P2"+"\t"+str(i+1)+"\t"+"PEI"+"\t"+pei_atype[i*11+j]+"\t"+str(i*11+j+1)+"\t0.000\n")
        else:
            atom_section.append(str(i*11+j+1)+"\t"+"Qd"+"\t"+str(i+1)+"\t"+"PEI"+"\t"+pei_atype[i*11+j]+"\t"+str(i*11+j+1)+"\t1.000\n")


bond_head="""
 
[ bonds ]
;  ai  aj   funct
"""


"""
for i in range(len(content)):
    if content[i][0:5]=="angle":
        print(i)
    elif content[i][0:5]=="dihed":
        print(i)
"""


info=open("Bond_length.txt")
content=info.readlines()
content[0].split()


bond={}
for i in range(1,13):
    #print(i,content[i].split(","))
    pair=content[i].split(",")[0]
    kb=content[i].split(",")[1]
    kd=content[i].split(",")[2]
    
    
    a,b=pair.split('-')
    bond[(a,b)]=[kd,kb]

angle={}
for i in range(14,49):
#    print(i,content[i].split(","))
    pair=content[i].split(",")[0]
    kb=content[i].split(",")[1]
    kd=content[i].split(",")[2]
    
    
    a,b,c=pair.split('-')
    angle[(a,b,c)]=[kd,kb]

angle_head="""
[ angles ]
;   ai     aj     ak    funct   theta         cth

"""

dihe_head="""
[ dihedrals ] ; propers
; for gromacs 4.5 or higher, using funct 9
;    i      j      k      l   func   phase     kd      pn
"""

#print(angle_head)
#print(dihe_head)

dihe_list=[]
for i in range(50,len(content)):
    if "s" in content[i] or "p" in content[i] or "t" in content[i]:
#        print(i)
        dihe_list.append(i)

dihe={}
for i in range(len(dihe_list)-1):
    a=content[dihe_list[i]].split("-")[0]
    b=content[dihe_list[i]].split("-")[1]
    c=content[dihe_list[i]].split("-")[2]
    d=content[dihe_list[i]].split("-")[3]
    d=d.split(',')[0]
    parm=[]
    for j in range(dihe_list[i]+1,dihe_list[i+1]):
        #print(j)
        phase=content[j].split(',')[0]
        kd=content[j].split(',')[1]
        pn=content[j].split(',')[2]
        parm.append(phase)
        parm.append(kd)
        parm.append(pn)
    dihe[(a,b,c,d)]=parm
    if ',' in d:
        print("....",d)


bond_section=[]
for i in range(len(conect)):
    a=pei_atype[conect[i][0]-1]
    b=pei_atype[conect[i][1]-1]
    
    if (a,b) in bond.keys():
        bond_section.append(str(conect[i][0])+"\t"+str(conect[i][1])+"\t1\t"+bond[(a,b)][0]+"\t"+bond[(a,b)][1]+"\n")
    elif (b,a) in bond.keys():
        bond_section.append(str(conect[i][0])+"\t"+str(conect[i][1])+"\t1\t"+bond[(b,a)][0]+"\t"+bond[(b,a)][1]+"\n")
    else:
        print("warning:",a,b)
        print(conect[i])


angle_section=[]
for i in range(len(anglelist)):
    a=pei_atype[anglelist[i][0]-1]
    b=pei_atype[anglelist[i][1]-1]
    c=pei_atype[anglelist[i][2]-1]
    
    if (a,b,c) in angle.keys():
        angle_section.append(str(anglelist[i][0])+"\t"+str(anglelist[i][1])+"\t"+str(anglelist[i][2])+"\t1\t"+angle[(a,b,c)][0]+"\t"+angle[(a,b,c)][1]+"\n")
    elif (c,b,a) in angle.keys():
        angle_section.append(str(anglelist[i][0])+"\t"+str(anglelist[i][1])+"\t"+str(anglelist[i][2])+"\t1\t"+angle[(c,b,a)][0]+"\t"+angle[(c,b,a)][1]+"\n")
    else:
        #print("warning:",a,b,c)
        print(anglelist[i])
        
dihe_section=[]
for i in range(len(dihelist)):
    a=pei_atype[dihelist[i][0]-1]
    b=pei_atype[dihelist[i][1]-1]
    c=pei_atype[dihelist[i][2]-1]
    d=pei_atype[dihelist[i][3]-1]
    
    if (a,b,c,d) in dihe.keys():
        x=dihe[(a,b,c,d)]
        for l in range((len(x)//3)):
            dihe_section.append(str(dihelist[i][0])+"\t"+str(dihelist[i][1])+"\t"+str(dihelist[i][2])+"\t"+str(dihelist[i][3])+"\t9\t"+x[l*3]+"\t"+x[l*3+1]+"\t"+x[l*3+2]+"\n")
    elif (d,c,b,a) in dihe.keys():
        x=dihe[(d,c,b,a)]
        for l in range((len(x)//3)):
            dihe_section.append(str(dihelist[i][0])+"\t"+str(dihelist[i][1])+"\t"+str(dihelist[i][2])+"\t"+str(dihelist[i][3])+"\t9\t"+x[l*3]+"\t"+x[l*3+1]+"\t"+x[l*3+2]+"\n")
    else:
        #print("warning:",a,b,c)
        print(dihelist[i])

output=open("pei.itp","w")
output.writelines(atom_head)

output.writelines(atom_section)

output.writelines(bond_head)

output.writelines(bond_section)

output.writelines(angle_head)

output.writelines(angle_section)

output.writelines(dihe_head)

output.writelines(dihe_section)

output.close()














