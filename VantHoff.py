__author__ = 'gvtheen'
import numpy as np
import copy
import ast
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from math import sqrt
class VantHoffFit(object):
    def __init__(self,filename):
        try:
            f=open(filename,'r')
        except IOError:
            print ('File not Exists')
        finally:
            f.close()
        self.filename=filename
        self.data=[]

    def readfile(self):
        fp=open(self.filename,'r')
        data=[0.0 for i in range(2) ]
        line=fp.readline()
        while(line!=''):
            line.strip()
            if line=='':
                break
            a=line.split()
            num=0
            for j in range(len(a)):
                if a[j]!="":
                    data[num]=ast.literal_eval(a[j])
                    num=num+1
                    if num==2:
                        break
            self.data.append(copy.deepcopy(data))
            line=fp.readline()

        self.data=np.array(self.data)
        self.data=self.data.transpose()
        #print(self.data[0,:])

    def fit(self):
        p0=[12500,-100,15.45, 2.45]
        XX=self.data[0,:]
        YY=2.828*np.sqrt(self.data[1,:])
        Para=leastsq(self.error,p0,args=(XX,YY))
        print("detH: %f" % Para[0][0])
        print("detS: %f" % Para[0][1])
        print("uffHt: %f" % Para[0][2])
        print("uffLt: %f" % Para[0][3])
        Te=Para[0][0]/Para[0][1]
        print("T1/2: %f" % Te)
        newX=self.data[0,:]
        newX=np.array(newX)

        newY=self.funcVantHoff(Para[0],newX)
        newFileName = "out_" + self.filename
        np.savetxt(newFileName,np.array([newX,(newY/2.828)*(newY/2.828)]).transpose(),fmt="%s")
        np.savetxt("para_"+self.filename,Para,fmt="%s")
        #np.savetxt("newY.txt",(newY/2.828)*(newY/2.828))
        plt.plot(XX,YY,newX,newY)
        plt.show()

    def funcVantHoff(self,p,T):
        detH,detS,uffHt,uffLt=p
        Tinv= 1/T
        R=8.314
        aa=np.exp(-1*detH*Tinv/R + detS/R)
        #uffHt=7.19
        #uffLt=4.75

        return np.sqrt((aa*uffHt*uffHt+uffLt*uffLt)/(aa+1))

    def error(self,p,x,y):
        detH,detS,uffHt,uffLt=p
        return self.funcVantHoff(p,x)-y

if __name__ == "__main__":
    aa = VantHoffFit("co3.txt")
    aa.readfile()
    aa.fit()

