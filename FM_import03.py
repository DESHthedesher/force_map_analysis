# -*- coding: utf-8 -*-
"""
Created on Sat May 17 13:24:46 2014

@author: Duncan
"""

#class which allows manipulation of workign directory
#and choice of what happens
import os
from matplotlib import pyplot as plt
from constants import *
import LS_Regress

##build things
LSR = LS_Regress.LS_Regress()      


class FM_UI():

    def __init__(self, *args, **kwargs):
        self.runner()

    def runner(self):
        
        curve_array = self.fmap_namer()
        
        for i in curve_array:
            defl = self.data_import(i[1])
            zsens = self.data_import(i[2])
            
            #checks that the curve isn't crap
            if defl.index(max(defl)) > 200:   
                defl, zsens = self.extend_readin(defl, zsens)
                defl = self.invols_correct(defl, zsens)
                force, sep = self.sep_and_force(defl, zsens)
        
        plt.plot(sep, force)
        plt.show()
                   

    #gets names of all complete ibw files in the folder
    #returns a 2D array containing all the filenames
    def fmap_namer(self):

        curve_array = []
        
        #gets all filenames in the folder containing Defl
        pathh = os.getcwd()
        listing = os.listdir(pathh)
        namelist = []
        for i in listing:
            if "Defl" in i:
                namelist.append(i)

    
        #iterates data_readin over all files in namelist
        for i in namelist:
            line = i[4:8]
            point = i[13:17]
            curvename = 'Line%sPoint%s' %(line,point)
            print '\n\n\n ***Reading Line%sPoint%s***\n' %(line,point)
            d = 'Line%sPoint%sDefl.txt' %(line,point)
            z = 'Line%sPoint%sZSnsr.txt' %(line,point)
            curve_ID = [curvename,d,z]
            curve_array.append(curve_ID)
        
        return curve_array


        
    def data_import(self,var):
        ###actually opens, reads, and close the nfile containing the data of interest
        if os.path.exists(var):
            fial = open(var)
            data = fial.read()
            data = data.split()
            data = [float(i) for i in data]
            fial.close()
            return data
        else:
            print 'Insufficient data for %s' %(var)
            data = ['no']
            return data
 

           
    def extend_readin(self, d, z):
        # makes lists of data of each type
        # separates the extension and retraction parts of the curve
        # should be upside down before invols_correct
        
        defl_ext = []
        zsen_ext = []

        #define where retraction begins as max in deflection + # data points/1000
        maxd = d.index(max(d))
        ret_point = maxd - len(d)/1000

        for i in range(0, ret_point):
            defl_ext.append(d[ret_point - i])
            zsen_ext.append(z[ret_point - i])
            
        return defl_ext, zsen_ext


        
    def retract_readin(self,d,z):
        # makes lists of data of each type
        # separates the extension and retraction parts of the curve
        
        defl_ret = []
        zsen_ret = []

        #define where retraction begins as max in deflection + # data points/1000
        maxd = d.index(max(d))
        ret_point = maxd + len(d)/1000

        for i in range(ret_point, len(d)):
            defl_ret.append(d[i])
            zsen_ret.append(z[i])

        return defl_ret, zsen_ret        
        

        
    def flat_average(self, d):
        #calculate an average deflection value using the 85 - 92 % of data points
        left_limit = int(LEFT_TAIL_LIMIT*len(d))
        right_limit = int(RIGHT_TAIL_LIMIT*len(d))
            
        average = 0.0
        for i in range(left_limit, right_limit):
            average = average + d[i]
        average = average/float(len(d)*0.17)
        
        return average
    


        
    def invols_correct(self, d, z):
        #converts values to DeflV using the exported invols, and re-calculations defl (nm)
        #based on the real value of invols as measured for each curve.

        #convert nm to V
        dV = []
        for i in d:
            new = i*INVOLS*10**5
            dV.append(new)

        #calculate an average deflV for the tail near the event
        avg = self.flat_average(d)

        #make a list of numbers from which invols may be calculated. uses top 1/3 of the curve
        count = dV.index(max(dV))
        indent_height = max(dV) - avg
        smallest_invols_point = indent_height/3.0

        involist = []
        while True:
            if dV[count] >= avg + smallest_invols_point:
                involist.append(count)
                count += 1
            else:
                break

        # find slope of contact region using least squares method
        #intercept, slope = LSR.linefit_noError(z,d,involist)
        
        #print "module",slope
        
        sumx = 0.0
        sumy = 0.00
        sumxx = 0.00
        sumxy = 0.0
        for i in involist:
            sumx += z[i]
            sumy += dV[i]
            sumxx += (z[i])**2
            sumxy += (z[i])*(dV[i])
        
        slope = (sumx*sumy - len(involist)*sumxy)/(sumx**2 - (len(involist))*sumxx)
        
        realinvols = (1/slope)*10**9

        # convert deflV back to nm using the correct value of realinVols
        d = []
        for i in dV:
            new = (i*float(realinvols))/(10**9)
            d.append(new)
         

        return d


    def sep_and_force(self,d,z):
        #converts deflection and zsensor to force and separation

        sep = []
        for i in range(0,len(z)):
            sep_val = -z[i] + d[i]
            sep.append(sep_val)        

        force = []
        for i in range(0, len(d)):
            f = SPRING_CONSTANT*(d[i])
            force.append(f)

        return force, sep
        

if __name__ == "__main__":
    FM_UI(None)
