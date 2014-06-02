# -*- coding: utf-8 -*-
"""
Created on Sat May 17 12:36:51 2014

@author: Duncan
"""
import math
from matplotlib import pyplot as plt

class indentation_fits():
    def __init__(self, *args, **kwargs):
        self.runner()
        
    ###tests all the functions    
    def runner(self):
        #a, b = self.JKR_gen(1,1,[1],0.5)
        #a = self.JKR_LS([1],[1,2,3,4,5,5],1,0.5)
        
        print "all good"

    def Exponent_gen(self, x, A, b, k):
        e = 2.71828
        
        fity = []
        for i in x:
            fity.append(A*e**(k*i) + b)

        return x, fity
        

    ###fits to the equation y = Ae^kx, assuming k is known
    def Exponent1_LS(self, x, y, k):
        e = 2.71828
        
        sumxy = 0.0
        sumxx = 0.0
        sumx = 0.0
        sumy = 0.0
        
        for i in range(0, len(x) - 1):
            sumx += e**(k*x[i])
            sumxy += y[i]*e**(k*x[i])
            sumxx += e**(2*k*x[i])
            sumy += y[i]
            
        b = (sumy*sumxx - sumxy*sumx)/((len(x)-1)*sumxx - sumx**2)
        A = (sumxy - b*sumx)/sumxx
        
        return A, b
        

    ###generates an idealized force distance profile based on the JKR contact
    ###model, using the parameters given
    def JKR_gen(self, E, R, x, a):
        
        fitx = []
        for i in range(0, len(x)):
            fitx.append(math.fabs(x[i] - x[-1]))
        
        force = []
        for i in fitx:
            new_force = (i**1.5)*((R**0.5)*E*4)/(3*(1-a**2))
            force.append(new_force)
        
        return x, force
    
    
    ### finds the best young's modulus for a given approach profile
    ### fits the data to this value and estimates a value for chi square (X2)
    ### x is distance, y is force, E is young's modulus, R is tip radius, a is poisson's ratio
    def JKR_LS(self,x,y,R,a):

        fitx = []
        for i in range(0, len(x)):
            fitx.append(math.fabs(x[i] - x[len(x)-1]))
            
        sumXY = 0.0
        sumXX = 0.0
        sumXandHalf = 0.0
        
        for i in range(0, len(x)):
            sumXX += fitx[i]**3
            sumXY += y[i]*(fitx[i]**1.5)
            sumXandHalf += fitx[i]**1.5
            
        A = 4*(R**0.5)/(3*(1-a**2))
        E = sumXY/(sumXX*A)
        
        #generate a fit based on calculated E
        fitforce = []
        for i in fitx:
            new_force = (i**1.5)*((R**0.5)*E*4)/(3*(1-a**2))
            fitforce.append(new_force)

        #calculate the average variance for the sample
        s2 = 0.0
        for i in range(0, len(y)):
            s2 += (y[i] - fitforce[i])**2
        
        s2 = s2/(len(y)- 2)

        #calculate the error in E
        Se = 0.0
        for i in fitx:
            Se += (sumXandHalf/sumXX)**2
        Se = ((Se*s2)/A)**0.5
       
        
        return E,Se
    
    ##generates a coulmbic fit profile
    def coulomb_gen(self, r, q):
        ke = 8.9875517873681764*10**9 ##coulomb's constant, N/m^2*C^2
        
        coulombic_force = []
        for i in range(0, len(r)):
            coulombic_force.append(ke*(q**2)/(r[i]**2))
        
        return r, coulombic_force
    
    ###finds the best |charge|(q) to fit the coulombic expression with
    ###assumes |charge| on each body is the same 
    def coulomb_LS(self,r,F):
        ke = 8.9875517873681764*10**9 ##coulomb's constant, N/m^2*C^2

        
        sumrf = 0.0
        sumrr = 0.0
        
        for i in range(0, len(F)):
            sumrf += (F[i])/(r[i]**2)
            sumrr += 1.0/(r[i]**4)
            
        q = ((1/ke)*(sumrf/sumrr))**0.5
        
        return q


    ###y is the data, VdW is from van der waals, EDL is from electric double layer
    def simplified_DLVO_LS(self, x, y, debye_length):
        e = 2.71828
        
        sumyVdW = 0.0
        sumVdWEDL = 0.0
        sumyEDL = 0.0
        sumEDLEDL = 0.0
        sumVdWVdW = 0.0
        
        for i in range(0, len(x) - 1):
            sumyVdW += y[i]*(1/(x[i]**2))
            sumVdWEDL += (1/(x[i]**2))*e**-(debye_length*x[i])
            sumyEDL += y[i]*e**-(debye_length*x[i])
            sumEDLEDL += (e**-(2*debye_length*x[i]))
            sumVdWVdW += (1/(x[i]**4))
            
        phi = (sumEDLEDL*sumyVdW - sumVdWEDL*sumyEDL)/(sumVdWEDL**2 - sumEDLEDL*sumVdWVdW)
        theta = (sumyEDL + phi*sumVdWEDL)/(sumEDLEDL)
        
        return phi, theta
        
        
    
    def simplified_DLVO_gen(self, debye_length, distance, theta, phi):
        e = 2.71828
        
        DLVO_forces = []
        for i in distance:
            nDLVO = theta*e**(-debye_length*i) - phi/(i**2)
            DLVO_forces.append(nDLVO)
            
        return DLVO_forces
        

if __name__ == "__main__":
    indentation_fits(None)