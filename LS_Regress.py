# -*- coding: utf-8 -*-
"""
Created on Thu May 22 16:01:47 2014

@author: Duncan
"""
class LS_Regress():
    
    def test(self):
        return "this is a test"

    #takes two arrays of floats and finds the least squares slope for points indexed at indexlist
    

    ###linear regression as outlined in bevvington 6.13
    ###fit to the equation y = a + bx
    ###takes x, y, and a list of the indexes for the points to be fitted
    ###returns intercept, slope
    
    def linefit_noError(self, x, y, indexlist):
        sumx = 0.0
        sumy = 0.0
        sumxy = 0.0
        sumxx = 0.0
        N = float(len(indexlist))

        for i in indexlist:
            sumx += x[i]
            sumxx += x[i]*x[i]
            sumy += y[i]
            sumxy += y[i]*x[i]
            
        delta = N*sumxx - sumx**2
        a = sumy*sumxx - sumx*sumxy
        b = sumxy*N - sumy*sumx
        
        a = a/delta
        b = b/delta
            
        return a, b

if __name__ == "__main__":
    LS_Regress(None)
