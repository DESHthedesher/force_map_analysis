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
import numpy as np
import math
import indentation_fits

##build things
LSR = LS_Regress.LS_Regress()      
indF = indentation_fits.indentation_fits()


class FM_UI():

    def __init__(self, *args, **kwargs):
        self.runner()

    def runner(self):
        
        curve_array = self.fmap_namer()
        
        for i in curve_array:
            print "\n\n\n****analyzing ",i[0],"*****"
            defl = self.data_import(i[1])
            zsens = self.data_import(i[2])
            
            #checks that the curve isn't crap
            if defl.index(max(defl)) > 200:   
                #processes data 
                defl, zsens = self.extend_readin(defl, zsens)
                defl = self.invols_correct(defl, zsens)
                force, sep = self.sep_and_force(defl, zsens)
                force, sep = self.baseline_correct(force, sep)
                force, sep = self.zero(force, sep)
                
                
                ##fits the region of interest to the model
                force_fit, sep_fit, indent = self.JKR_fitmap(force, sep, zsens, FIT_THE_DROP)
                reduced_force = []
                reduced_sep = []
                for i in range(0, len(force_fit)):
                    reduced_force.append(force_fit[i])
                    reduced_sep.append(sep_fit[i])
                    
                
                
                if indent == True:
                    print "indentation"
                elif indent != True:
                    print "no indentation"
                
                if MIN_FIT_LENGTH < len(force_fit) < 20*SLIDING_WINDOW_SIZE:
                    ##sneddon model
                    youngs_modulus, Eerr = indF.JKR_LS(reduced_sep, reduced_force, TIP_RADIUS, POISSON_RATIO)
                    JKR_sep, JKR_force = indF.JKR_gen(youngs_modulus, TIP_RADIUS, reduced_sep ,POISSON_RATIO)
                    
                    ##simple exponential
                    A, b = indF.Exponent1_LS(reduced_sep, reduced_force, -1*DEBYE_LENGTH)
                    exp_sep, exp_force = indF.Exponent_gen(reduced_sep,A,b, -1*DEBYE_LENGTH)

                    ##coulombic
                    q = indF.coulomb_LS(sep_fit, force_fit)
                    cool_sep, cool_force = indF.coulomb_gen(sep_fit, q)
                    
                    ##DLVO
                    phi, theta = indF.simplified_DLVO_LS(sep_fit, force_fit, DEBYE_LENGTH)
                    DLVO_force = indF.simplified_DLVO_gen(DEBYE_LENGTH, sep_fit, theta, phi)
                    c1, c2, c3 = indF.Nonzeroed_DLVO_LS(sep_fit, force_fit, DEBYE_LENGTH)
                    DLVO2_force = indF.Nonzeroed_DLVO_gen(sep_fit, c1,c2,c3,DEBYE_LENGTH)
                    print "A =", c2*6.0/TIP_RADIUS
                    print "Z =", c3/(TIP_RADIUS*DEBYE_LENGTH)
                    print "1/k =", 1/DEBYE_LENGTH
                    
                    print "E =",youngs_modulus," +/- =",Eerr
                    plt.plot(sep, force)
                    plt.plot(sep_fit,force_fit)
                    #print "q =",q
                    #plt.plot(cool_sep, cool_force)
                    #plt.plot(JKR_sep, JKR_force)
                    #plt.plot(exp_sep, exp_force)
                    #plt.plot(sep_fit, DLVO_force)
                    plt.plot(sep_fit, DLVO2_force)
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
        

        
    def flat_stats(self, y):
        #calculate an average deflection value using the 85 - 92 % of data points
        left_limit = int(LEFT_TAIL_LIMIT*len(y))
        right_limit = int(RIGHT_TAIL_LIMIT*len(y))
            
        average = 0.0
        for i in range(left_limit, right_limit):
            average = average + y[i]
        average = average/float(len(y)*0.17)
        
        variance = 0.0
        for i in range(left_limit, right_limit):
            variance += (average - y[i])**2
        variance = variance/(right_limit - left_limit)
        
        return average, variance
    


        
    def invols_correct(self, d, z):
        #converts values to DeflV using the exported invols, and re-calculations defl (nm)
        #based on the real value of invols as measured for each curve.

        #convert nm to V
        dV = []
        for i in d:
            new = i*INVOLS*10**5
            dV.append(new)

        #calculate an average deflV for the tail near the event
        avg, variance = self.flat_stats(d)

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
        intercept, slope = LSR.linefit_noError(z,d,involist)
        #realinvols = slope*10**9
        
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
    
    
    
    def find_nearest(self,array,value):
        
        nearest = max(array)
        for i in array:
            if math.fabs(i - value) < nearest:
                nearest = i
        
        return array.index(nearest)
    
    
    
    def zero(self, force, sep):
        ##places the curve such that Defl away from the surface is zero
        ##       and sep at the surface is zero
        
        ##zero the force
        forceaverage, variance = self.flat_stats(force)
        
        newforce = []
        for i in force:
            newforce.append(i - forceaverage) 
        
        
        ##zero the sep
        smallest_point = self.find_nearest(newforce, max(newforce)/3.0)
        
        total = 0
        for i in range(0, smallest_point):
            total += sep[i]
 
        sepaverage = total/(smallest_point)
        
        newsep = []
        for i in sep:
            newsep.append(i - sepaverage)
            
        return newforce, newsep
        
    def three_plots(self, y1, x1):
        ##Normalizes and plots three things at once
        ##max[y] - min[y] = 1.0
        print "Plotting...."
        
        y2, x2 = self.differentiation(y1, x1, SLIDING_WINDOW_SIZE)
        y3, x3 = self.differentiation(y2, x2, SLIDING_WINDOW_SIZE)
        
        
        a1 = 1.0/(max(y1)-min(y1))
        a2 = 1.0/(max(y2)-min(y2))
        a3 = 1.0/(max(y3)-min(y3))
        
        norm_y1 = []
        for i in y1:
            norm_y1.append(i*a1)
            
        norm_y2 = []
        for i in y2:
            norm_y2.append(math.fabs(i*a2 + 0.2))
            
        norm_y3 = []
        for i in y3:
            norm_y3.append(math.fabs(i*a3 + 1.3))
            
        norm_x1 = []
        for i in x1:
            norm_x1.append(-1*i)
        
        norm_x2 = []
        for i in x2:
            norm_x2.append(-1*i)
            
        norm_x3 = []
        for i in x3:
            norm_x3.append(-1*i)
 
        
        plt.plot(norm_x1, norm_y1)
        plt.plot(norm_x2, norm_y2)
        plt.plot(norm_x3, norm_y3)
        plt.show()
        
    def baseline_correct(self, force, sep):
        #Corrects for slope in the baseline.
        # Fcorrect  =  F initial - S*slope
        left_limit = int(LEFT_TAIL_LIMIT*len(force))
        right_limit = int(RIGHT_TAIL_LIMIT*len(force))
        
        pointlist = []
        for i in range(left_limit, right_limit):
            pointlist.append(i)
        
        intercept, slope = LSR.linefit_noError(sep, force, pointlist)
        
        for i in range(0, len(force)):
            newforce = force[i] - slope*sep[i]
            force[i] = newforce
        
        return force, sep
    
    
    def differentiation(self, y, x, window_size):
        #creates a second curve of df/ds for the entire curve
        #this is calculated as the average change to the left and right of
        #a given data point. The average is taken over the value of SLIDING_WINDOW_SIZE
        left_limit = int(LEFT_TAIL_LIMIT*len(x))
        
        dydx = []
        smaller_x = []
        for i in range(window_size, left_limit):
            
            left_average_y = 0.0
            left_average_x = 0.0
            for k in range(i - window_size, i):
                left_average_y += y[k]
                left_average_x += x[k]
            left_average_y = left_average_y/window_size
            left_average_x = left_average_x/window_size
            
            right_average_y = 0.0
            right_average_x = 0.0
            for k in range(i, i + window_size):
                right_average_y += y[k]
                right_average_x += x[k]
            right_average_y = right_average_y/window_size
            right_average_x = right_average_x/window_size
            
            rise_over_run = (right_average_y - left_average_y)/(right_average_x - left_average_x)
            dydx.append(rise_over_run)
            smaller_x.append(x[i])
            
        return dydx, smaller_x
        
        
    def JKR_fitmap(self, force, sep, zsens, fitregion):
        # finds region of interest for JKR fit.
        # defines region of interest as the region from where |slope| > 2stdev
        # to where dF/ds changes sign
        # fit the drop indicates that you want to fit the region where force begins to drop
        Indent_Detect = True
        
        dfdz, small_zsens = self.differentiation(force, zsens, SLIDING_WINDOW_SIZE)
        d2fdz2, small_zsens2 = self.differentiation(dfdz, small_zsens, SLIDING_WINDOW_SIZE)
        average, variance = self.flat_stats(dfdz)
        d1_average, d1_variance = self.flat_stats(dfdz)
        
        ##find where the tips starts to interact with the sample
        right_limit = len(dfdz) - 1
        number_of_outliers = 0
        while number_of_outliers < DETECTION_LIMIT:
            right_limit += -1
            
            if right_limit < 0:
                number_of_outliers = DETECTION_LIMIT
                right_limit = 1
            elif dfdz[right_limit] > average + (variance**0.5)*DETECTION_LIMIT:
                number_of_outliers += 1
            else:
                number_of_outliers = 0
        
        ##find where a breakthrough happens
        #left_limit = d2fdz2.index(min(d2fdz2))
        left_limit = max(d2fdz2)
        for i in range(1, min(right_limit, len(d2fdz2) - 1)):
            if math.fabs(d2fdz2[i]-d2fdz2[i-1]) > left_limit:
                left_limit = d2fdz2[i]
        left_limit = d2fdz2.index(left_limit)
        
        
        ##add some borders to the fits so that they are in region of interest only        
        if fitregion == True:
            leftFitIndex = SLIDING_WINDOW_SIZE*2 + left_limit + 1
        elif fitregion == False:
            leftFitIndex = SLIDING_WINDOW_SIZE*3 + left_limit + 1
        rightFitIndex = SLIDING_WINDOW_SIZE*2 + right_limit + 1
        
        
        ## if the fit is too short, try using first derivative maximum to find it (no breakthrough)
        if rightFitIndex - leftFitIndex < MIN_FIT_LENGTH:
            Indent_Detect = False
            left_limit = max(dfdz)
            tempLeftLimit = max(dfdz)
            
            for i in range(1, min(right_limit, len(dfdz) - 1)):
                if left_limit - dfdz[i] > left_limit - DETECTION_LIMIT*d1_variance**0.5:
                    tempLeftLimit = dfdz[i]
             
            
            leftFitIndex = dfdz.index(tempLeftLimit)

        
        force_fit = []
        sep_fit = []
        for i in range(leftFitIndex, rightFitIndex):
            force_fit.append(force[i])
            sep_fit.append(sep[i])
            
        return force_fit, sep_fit, Indent_Detect
        
        
 
    
            

if __name__ == "__main__":
    FM_UI(None)
