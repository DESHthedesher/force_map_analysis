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
    
    
    def fit_variance(self, realY, fitY, Nparams):
        
        s2 = 0.0
        for i in range(0, len(realY) - 1):
            s2 += (realY[i] - fitY[i])**2
        
        s2 = s2/(len(realY)- 1 - Nparams)
    
        return s2
    
    #calculates Chi square for data with a given fit and variance
    def ChiSquare(self, realy, fity, variance):
        
        X2 = 0.0
        for i in range(0, len(variance) - 1):
            X2 += ((realy[i]-fity[i])**2)/variance[i]
        
        return X2
       
    #calculates Chi Square for data with a given fit and unknown variance   
    def ChiSquare_UnkownErr(self, realy, fity, Nparams):
        
        v = self.fit_variance(realy, fity, Nparams)
        variance = []
        for i in range(0, len(realy) - 1):
            variance.append(v)
            
        X2 = self.ChiSquare(realy, fity, variance)
        return X2
    
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
