# -*- coding: utf-8 -*-
"""
Created on Sat May 17 12:36:51 2014

@author: Duncan
"""
import numpy as np
from matplotlib import pyplot as plt

class JKR():
    def __init__(self, *args, **kwargs):
        self.runner()
        
    ###tests all the functions    
    def runner(self):
        a, b = self.JKR_gen(1,1,[1],0.5)
        a = self.JKR_LS([1],[1],1,0.5)
        
        print "all good"

    ###generates an idealized force distance profile based on the JKR contact
    ###model, using the parameters given
    def JKR_gen(self, E, R, x, a):
        for i in range(0, len(x)):
            if x[i] < 0:
                x[i] = 0
        
        force = []
        for i in x:
            force.append(((4*E*(R**0.5))/(3*(1-a**2)))*(i**1.5))
        
        return x, force
    
    
    ### finds the best young's modulus for a given approach profile
    ### x is distance, y is force, E is young's modulus, R is tip radius, a is poisson's ratio
    def JKR_LS(self,x,y,R,a):

        for i in range(0, len(x)):
            if x[i] < 0:
                x[i] = 0
        print x
        
        sumXY = 0.0
        sumXX = 0.0
        
        for i in range(0, len(x)):
            sumXX += x[i]**3
            sumXY += y[i]*(x[i]**1.5)
            
        E = ((3*(1-a**2))/(4*R**0.5))*(sumXY/sumXX)
        
        return E
        

if __name__ == "__main__":
    JKR(None)