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

    ###generates an idealized force distance profile based on the JKR contact
    ###model, using the parameters given
    def JKR_gen(self, E, R, x, a):
        
        fitx = []
        for i in range(0, len(x)):
            fitx.append(math.fabs(x[i] - x[-1]))
        
        force = []
        for i in fitx:
            #force.append(1000*((4*E*(R**0.5))/(3*(1-a**2)))*(i**1.5))
            new_force = (i**1.5)*((4.0/3.0)*E*R**0.5)/(3-3*a**2)
            force.append(new_force)
        
        return x, force
    
    
    ### finds the best young's modulus for a given approach profile
    ### x is distance, y is force, E is young's modulus, R is tip radius, a is poisson's ratio
    def JKR_LS(self,x,y,R,a):

        newx = []
        for i in range(0, len(x)):
            newx.append(math.fabs(x[i] - x[-1]))
            
        sumXY = 0.0
        sumXX = 0.0
        
        for i in range(0, len(x)):
            sumXX += newx[i]**3
            sumXY += y[i]*(newx[i]**1.5)
            
        nom = ((3*(1-a**2))/(4*R**0.5))
        denom = (sumXY/sumXX)
        E = nom/denom
        return E
    
    ##generates a coulmbic fit profile
    def coulomb_gen(self, r, q):
        ke = 8.987*10**9 ##coulomb's constant, N/m^2*C^2
        
        coulombic_force = []
        for i in range(0, len(r)):
            coulombic_force.append(ke*(q**2)/(r[i]**2))
        
        return r, coulombic_force
    
    ###finds the best |charge|(q) to fit the coulombic expression with
    ###assumes |charge| on each body is the same 
    def coulomb_LS(self,r,F):
        ke = 8.987*10**9 ##coulomb's constant, N/m^2*C^2
        
        sumrf = 0.0
        sumrr = 0.0
        
        for i in range(0, len(F)):
            sumrf += (F[i]**2)/(r[i]**2)
            sumrr += 1.0/(r[i]**4)
            
        q = (ke*sumrf/sumrr)**0.5
        
        return q
        
        

if __name__ == "__main__":
    indentation_fits(None)