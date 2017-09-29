import numpy
from numpy import linalg 
import scipy.linalg
import random
import math
#Let's define the molecule class that should contain the information in the format below

#5
#0001 -417.031
#C      1.04168000 -0.05620000 -0.07148000    1.04168200 -0.05620000 -0.07148100
#H      2.15109000 -0.05620000 -0.07150000    2.13089400 -0.05620200 -0.07149600
#H      0.67187000  0.17923000 -1.09059000    0.67859800  0.17494100 -1.07204400
#H      0.67188000  0.70866000  0.64196000    0.67861300  0.69474600  0.62898000
#H      0.67188000 -1.05649000  0.23421000    0.67861400 -1.03828500  0.22864100

class Molecule:
    def __init__(self, id, nat, atomization):
        self.id = id
        self.nat = nat
        self.atomization = atomization
        self.nh = 0 
        self.logic = 0
        self.coulomb = []
        self.atoms = []

    def setSelected(self, value):
        self._selected = value
    

def charge(element):
    if (element=='C'):
        return 6
    elif (element=='O'): 
        return 8
    elif (element=='H'):
        return 1
    elif (element=='N'):
        return 7
    elif (element=='S'):
        return 16
    else:
        print "element not identified", element

def number_of_h(nat, typ):

    loopa=0
    n_h=0
    while (loopa < nat):
        if (typ[loopa]=='H'):
            n_h=n_h+1
        loopa=loopa+1
    return n_h

def coulmatrix(nat, x, y, z, typ):
   
#    print "calling me"
    loopx=0
    cm=numpy.zeros((23,23))
#    print nat, x, y, z, typ
    while (loopx < nat):

        zi=charge(typ[loopx])
        cm[loopx][loopx]=0.5*zi**(2.4)
        loopy=0        

        while (loopy < nat):


            if (loopx != loopy): 
                zj=charge(typ[loopy])
                distance=( (float(x[loopx])-float(x[loopy]))**2 + (float(y[loopx])-float(y[loopy]))**2 + (float(z[loopx])-float(z[loopy]))**2  )**0.5
                cm[loopx][loopy]=float(zi)*float(zj)/distance
            loopy=loopy+1

       
        loopx=loopx+1

#    print "cm",cm
    normarray=[]
    loopx=0
    while (loopx < nat):

        tmpnorm=numpy.linalg.norm(cm[loopx][:])
        normarray.append(tmpnorm)
        loopx=loopx+1

    loopx=0
    posarray=[]
    while (loopx < nat):
        loopy=0
        position=0
        while (loopy < nat):
            if (float(normarray[loopx])<float(normarray[loopy])):
                position=position+1
            elif (float(normarray[loopx])==float(normarray[loopy]) and loopx>loopy):
                position=posarray[loopy]+1
                loopy=1000 
            loopy=loopy+1
        loopx=loopx+1
        posarray.append(position)


    loopx=0
    cm2=numpy.zeros((23,23))
#    print nat, x, y, z, typ
    while (loopx < nat):
        loopy=0
        loopx_tmp=posarray[loopx]

        while (loopy < nat):

            loopy_tmp=posarray[loopy]
#            print "loop",loopx_tmp, loopy_tmp
            cm2[loopx_tmp][loopy_tmp]=cm[loopx][loopy]
            loopy=loopy+1

        loopx=loopx+1

#    print "cm2",cm2

#    print "attenzione0"
#    print normarray
#    print posarray
#    print "attenzione1"
    loopx=0
    loopz=0
    v=numpy.zeros((276))
#    print nat, x, y, z, typ
    while (loopx < 23):
        loopy=loopx

        while (loopy < 23):

            v[loopz]=cm2[loopx][loopy]
            loopy=loopy+1
            loopz=loopz+1

        loopx=loopx+1        

#    print "loopz",loopz
#    print v

#    print cm
#    v, w=linalg.eigh(cm)
#    print v
#    v_ord=sorted(v, reverse=True, key=abs)
#    v_to_sort=abs(v)
#    v_ord=sorted(v_to_sort, reverse=True)
#    v_out=numpy.asarray(v_ord)
#    print v_out
    return v

moldata = open ('dsgdb7ae2.xyz','r')
infile = moldata
moldata.close

classarray = []
loop0=0
nat_smallh=0
while loop0 < 7102:
 
    line1 = infile.readline()
    line2 = infile.readline()
    tmpsplit0=line2.split()
    clmol=Molecule(loop0,line1,tmpsplit0[1])
#    print line1
    loop=0
    x=[]
    y=[]
    z=[]

    while loop < int(line1):
        line3 = infile.readline()
        tmpsplit=line3.split()
        clmol.atoms.append(tmpsplit[0])
#        print charge(tmpsplit[0])
        x.append(tmpsplit[1])
        y.append(tmpsplit[2])
        z.append(tmpsplit[3])
        loop=loop+1

    vv = coulmatrix(int(line1),x,y,z,clmol.atoms)
    nhh = number_of_h(int(line1),clmol.atoms)
    clmol.coulomb=vv
    clmol.nh=nhh
#    print clmol.nh
    if ((int(line1)-nhh)<5): 
        clmol.logic=1
        nat_smallh=nat_smallh+1
#    print "lllll",clmol.nh

        
#        clmol.x.append(tmpsplit[1])
#        clmol.y.append(tmpsplit[2])
#        clmol.z.append(tmpsplit[3])
#        clmol.xdf.append(tmpsplit[4])
#        clmol.ydf.append(tmpsplit[5])
#        clmol.zdf.append(tmpsplit[6])
    
    classarray.append(clmol) 
    loop0=loop0+1

kernel=numpy.zeros((1000,1000))

#print nat_smallh

loop0=0
while loop0 < (900-nat_smallh):
    index = random.randint(0,7101)
    if (classarray[index].logic == 0):
        classarray[index].logic=1
        loop0=loop0+1
#        print index

loop0=0
while loop0 < 100:
    index = random.randint(0,7101)
    if (classarray[index].logic == 0):
        classarray[index].logic=2
        loop0=loop0+1
#        print index

#    clmol.nh=nhh    

kernel=numpy.zeros((900,900))
yyy=numpy.zeros((900))

loop0=0
m1=0
while loop0 < 7102:
    if (classarray[loop0].logic==1):
        m2=0
        loop1=0
        yyy[m1]=classarray[loop0].atomization
#        print classarray[loop0].coulomb
        while loop1 < 7102:
            if (classarray[loop1].logic==1):
                distance=numpy.linalg.norm(classarray[loop0].coulomb-classarray[loop1].coulomb,ord=1)
                kernel[m1][m2]=math.exp( -distance/(4096) )
                if (m1==m2):
                    kernel[m1][m2]=kernel[m1][m2]+0.00000000093132257
                m2=m2+1
            loop1=loop1+1
        m1=m1+1
    loop0=loop0+1
#        ksksk classarray


cho = scipy.linalg.cho_factor(kernel)
alpha = scipy.linalg.cho_solve(cho, yyy)
#print alpha

loop0=0
m1=0
elle=numpy.zeros((900,100))
while loop0 < 7102:
    if (classarray[loop0].logic==1):
        m2=0
        loop1=0
#        print classarray[loop0].coulomb
        while loop1 < 7102:
            if (classarray[loop1].logic==2):
                distance=numpy.linalg.norm(classarray[loop0].coulomb-classarray[loop1].coulomb,ord=1)
                elle[m1][m2]=math.exp( -distance/(4096) )
#                if (m1==m2):
#                    elle[m1][m2]=elle[m1][m2]+0.000000000001
                m2=m2+1
            loop1=loop1+1
        m1=m1+1
    loop0=loop0+1

prediction=numpy.dot(elle.T,alpha)
#print prediction

loop1=0
m1=0
mae=0.0
while loop1 < 7102:
            if (classarray[loop1].logic==2):
                mae=mae+math.fabs(float(classarray[loop1].atomization)-prediction[m1])
                m1=m1+1
            loop1=loop1+1

print "mae",mae/100
#print clmol.y
#print clmol.zdf
#print classarray[0].zdf, classarray[1].zdf, classarray[99].zdf

loop0=0
m1=0
elle=numpy.zeros((900,6102))
while loop0 < 7102:
    if (classarray[loop0].logic==1):
        m2=0
        loop1=0
#        print classarray[loop0].coulomb
        while loop1 < 7102:
            if (classarray[loop1].logic==0):
                distance=numpy.linalg.norm(classarray[loop0].coulomb-classarray[loop1].coulomb,ord=1)
                elle[m1][m2]=math.exp( -distance/(4096) )
#                if (m1==m2):
#                    elle[m1][m2]=elle[m1][m2]+0.000000000001
                m2=m2+1
            loop1=loop1+1
        m1=m1+1
    loop0=loop0+1

prediction=numpy.dot(elle.T,alpha)
#print prediction

loop1=0
m1=0
mae=0.0
while loop1 < 7102:
            if (classarray[loop1].logic==0):
                mae=mae+math.fabs(float(classarray[loop1].atomization)-prediction[m1])
                m1=m1+1
            loop1=loop1+1

print "final mae",mae/6102
