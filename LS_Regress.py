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
        N = len(indexlist)

        for i in indexlist:
            sumx += x[i]
            sumy += y[i]
            sumxy += x[i]*y[i]
            sumxx += (x[i])**2

        delta = N*sumxx - sumx**2
        a = (1/delta)*(sumy*sumxx - sumx*sumxy)
        b = (1/delta)*(N*sumxx - sumx*sumy)

	b = (sumx*sumy - sumxx*N)/(sumx**2 - sumxx*N)

	return a, b

if __name__ == "__main__":
    LS_Regress(None)
