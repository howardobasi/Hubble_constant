#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:44:33 2020

@author: howardobasi
"""

import numpy as np 
from matplotlib import pyplot as plt 
import scipy.optimize as opt
import math as m 

def cpfunc(a,b,c):
    '''param: a 'dummy' gradient 
              b  changing variable in the function
              c  'dummy' y intercept 
       return: dummy function for opt.curvefit '''
    return a * b + c

def chi2(y,ysig,ym):
    '''param: y - observed y values 
              ysig - error on the observed y value 
              ym - expected y values 
       dtype: (array of) floats 
       return: chi squared value of the data'''
    chi2 = np.sum(np.square((y-ym)/(ysig)))
    return chi2

def corelation(x,y,z):
    '''param: x - covariance between the parameters 
              y - absolute error on the gradient 
              z - absolute error on the y intercept 
       dtype: float
       return: correlation coefficient of the errors on the two parameters'''
    a=x/(y*z)
    return a

def gerror(y):
    '''param: y - variance of the parameter      
       dtype: float
       return: standard deviation of the parameter '''
    x=np.sqrt(y)
    return x

def distfunc (d,x,y,z,c,a,b):
 '''param: d - log(period) of the galaxy
           x - apparent magnitude of the galaxy
           y - alpha constant (gradient parameter)
           z - beta constant (y intercept parameter)
           c - extinction of the galaxy
           a - alpha const error 
           b - beta const error 
    dtype : float (array)
    return : gal_dist_avg - average distance from earth of the galaxy 
             dist_err_avg - average distance from earth error of the galaxy
             gal_dist_err - array of the distances from earth for all stars in said galaxy
             gal_dm_err -  array of the distance from earth errors for all stars in said galaxy
          
'''
 gal_ab = (y * d) + z
 gal_dm = x - gal_ab
 gal_dm2 = (gal_dm + 5 - c)/5
 gal_dist = 10**(gal_dm2)
 gal_dist_avg = np.sum(gal_dist)/len(d)
 
 gal_ab_err = np.sqrt(np.square((a * d)) + np.square(b))     
 gal_dm_err = gal_ab_err      
 gal_dm_err2 = gal_dm_err/5
 gal_dist_err = gal_dist * np.log(10)* gal_dm_err2     
 dist_err_avg = np.sum(gal_dist_err)/len(d)

 
 return gal_dist_avg,dist_err_avg,gal_dist_err,gal_dm_err


def errorfit(x,y,z,a):
 '''param: x - array of the distance from earth errors for all stars in said galaxy
    dtype: float
    return: histogram of the errors to show distribution '''
 plt.xlim(x)
 plt.hist(y,bins=z,color="blue")
 a
 plt.ylabel('Number of Points')   
 plt.title('Error distribution ')
 
 return plt.show()

def gfunc(a,x):
    '''param: a - 'dummy' variable for the gradient parameter'
              x - 'dummy' variable for the changing variable in the function
       return: return: dummy function for opt.curvefit '''
    return  x * a

print('Assignment 2 - Hubbles Constant - Relevant Values')
print()
    
plax,plax_err,period ,ap_mag ,extin, extin_err = np.loadtxt('MW_Cepheids.dat',\
unpack=True,usecols=(1,2,3,4,5,6),dtype=float)  #inputing files with the relevant data and assinging variables to each column  

dfe = 1000/plax # equation relating the distance of the star from earth and the parralax 

dfe_err = (plax_err/(plax*plax))*1000 # propagating parralax error to the error on distance distance of the star from earth 

d_mod = (5 * np.log10(dfe)) - 5 + extin # equation relating distance modulus, distance and extinction

var1_err = 5 * (dfe_err/(dfe*np.log(10))) 
d_mod_err = np.sqrt(np.square(var1_err) + np.square(extin_err))  # propagating distance error and extinction error distance modulus error 

ab_mag = ap_mag - d_mod # euqation relating absoloute magnitude, apparent magnitude and distance modulus 

ab_mag_err = d_mod_err # error stays constant through above equation 

per_d = np.log10(period) # logging the period to get it into correct form of the equatiom

param, pcov = opt.curve_fit(f=cpfunc, xdata = per_d, ydata = ab_mag,\
sigma = d_mod_err,absolute_sigma=True) # finding the optimal parameters alpha and beta and their corresponding errors for the data 
  
print('parameters') 
print(param)
grad = param[0]         # printing the corresponding parameters computed from opt.curve_fit - alpha and beta 
inter = param[1]

y_new = grad * per_d + inter  # finding the expected y values using the optimal parameters, this will be used to compute chi^2

plt.ylim(-5,0)
plt.xlim(0.5,1.2)
plt.xlabel('log(P)')          #  plotting the graph of the log of the period against the absolute magnitude 
plt.ylabel('M')           
plt.title('Cepheid Period-Luminosity Relation')
plt.plot(per_d,ab_mag , color='red', marker='o', markersize=6, linestyle="None")
plt.errorbar(per_d, ab_mag, yerr=ab_mag_err, color='black', linestyle="None")
plt.plot(per_d, y_new, color='blue', marker="None", linewidth=2, linestyle="-") # plotting the trend line using the optimal parameters computed through opt.curvefit  
plt.show()

print()

grad_err = gerror(pcov[0,0])  #  computing the error on the two parameters 
inter_err = gerror(pcov[1,1])

print('a = ' ,grad,'+/-',grad_err)          #  printing the two parameters to show their values 
print('b = ' ,inter,'+/-',inter_err)
print()

cor_1 = corelation(pcov[0,1],grad_err,inter_err)  #  finding if the errors on the two parameters are correlated 
cor_2 = corelation(pcov[1,0],grad_err,inter_err)  

print('correlation coefficient =',cor_1,',',cor_2)  # correlation coefficent of errors on the two parameters 

print()
print('X^2 =',chi2(ab_mag,ab_mag_err,y_new))
print('Reduced X^2 =',chi2(ab_mag,ab_mag_err,y_new)/(len(per_d) -2))    # computing chi^2 and the reduced chi^2 

val = m.sqrt(2/(len(plax)-2))
print()
print('We expect reduced x^2 ~ 1 +/-',val) # computing what value we expect reduecd chi^2 to be around for a good fit 

mean = (np.sum(per_d))/len(per_d) 
new_per_d = per_d - mean    # as the two errors of the parameters are strongly correlated i am uncorrelating them by setting mean of the x values to 0 
print()
print('Taking log(P) mean = 0 ~', (np.sum(new_per_d))/len(new_per_d) ) 

print()

param2, pcov2 = opt.curve_fit(f=cpfunc, xdata = new_per_d, ydata = ab_mag,\
sigma=ab_mag_err,absolute_sigma=True) # finding the new optimal parameters alpha and beta and their corresponding errors for the data 

print('parameters')
print(param2)
grad2 = param2[0]       # printing the new corresponding parameters computed from opt.curve_fit - alpha and beta
inter2 = param2[1]
    
print()

y_2new=grad2 * new_per_d + inter2  # finding the expected y values using the optimal parameters, this will be used to compute chi^2

plt.ylim(-5,0)
plt.xlim(-0.5,0.6)          
plt.xlabel('log(P)')          #  plotting the graph of the new log of the period against the absolute magnitude 
plt.ylabel('M')
plt.title('Cepheid Period-Luminosity Relation (uncorrelated error)')
plt.plot(new_per_d,ab_mag , color='red', marker='o', markersize=6, linestyle="None")
plt.errorbar(new_per_d, ab_mag, yerr=ab_mag_err, color='black', linestyle="None")
plt.plot(new_per_d, y_2new, color='blue', marker="None", linewidth=2, linestyle="-")
plt.show()

grad_err2 = gerror(pcov2[0,0])    #  computing the error on the two new parameters
inter_err2 = gerror(pcov2[1,1])


print('a = ' ,grad2,'+/-',grad_err2)       #  printing the two new parameters to show their values
print('b = ' ,inter2,'+/-',inter_err2)
print()

cor2_1 = corelation(pcov2[0,1],grad_err2,inter_err2)
co2r_2 = corelation(pcov2[1,0],grad_err2,inter_err2)     # finding if the errors on the two new parameters are correlated this will tell us if my method worked to uncorrelate the two errors 

print('correlation coefficient =',cor2_1,',',co2r_2)


print()
print('X^2 =',chi2(ab_mag,ab_mag_err,y_2new))              # computing new chi^2 and the reduced chi^2 
print('Reduced X^2 =',chi2(ab_mag,ab_mag_err,y_2new)/8)

print()

gal_name = np.loadtxt('galaxy_data.dat',unpack=True,usecols=(0),dtype=str)  
v_rec,extin = np.loadtxt('galaxy_data.dat',unpack=True,usecols=(1,2),dtype=float)   # inputing relevant data and assigning a variable name to the columns 


gal1_p, gal1_ap =np.loadtxt('hst_gal1_cepheids.dat',unpack=True,usecols=(1,2),dtype=float)
gal2_p, gal2_ap =np.loadtxt('hst_gal2_cepheids.dat',unpack=True,usecols=(1,2),dtype=float)
gal3_p, gal3_ap =np.loadtxt('hst_gal3_cepheids.dat',unpack=True,usecols=(1,2),dtype=float)
gal4_p, gal4_ap =np.loadtxt('hst_gal4_cepheids.dat',unpack=True,usecols=(1,2),dtype=float)    # inputing more data and assinging variable names to each column 
gal5_p, gal5_ap =np.loadtxt('hst_gal5_cepheids.dat',unpack=True,usecols=(1,2),dtype=float)
gal6_p, gal6_ap =np.loadtxt('hst_gal6_cepheids.dat',unpack=True,usecols=(1,2),dtype=float)
gal7_p, gal7_ap =np.loadtxt('hst_gal7_cepheids.dat',unpack=True,usecols=(1,2),dtype=float)
gal8_p, gal8_ap =np.loadtxt('hst_gal8_cepheids.dat',unpack=True,usecols=(1,2),dtype=float)
    
gal1_dist,gal1_de,gal1_dist_err,gal1_dm_err = distfunc(gal1_p,gal1_ap,grad2,inter2,extin[0],grad_err2,inter_err2) 
gal2_dist,gal2_de,gal2_dist_err,gal2_dm_err = distfunc(gal2_p,gal2_ap,grad2,inter2,extin[1],grad_err2,inter_err2)
gal3_dist,gal3_de,gal3_dist_err,gal3_dm_err = distfunc(gal3_p,gal3_ap,grad2,inter2,extin[2],grad_err2,inter_err2)
gal4_dist,gal4_de,gal4_dist_err,gal4_dm_err = distfunc(gal4_p,gal4_ap,grad2,inter2,extin[3],grad_err2,inter_err2)   # using the function distfunc to work out distance for each galaxy and assignign a variable name to eahc output
gal5_dist,gal5_de,gal5_dist_err,gal5_dm_err = distfunc(gal5_p,gal5_ap,grad2,inter2,extin[4],grad_err2,inter_err2)
gal6_dist,gal6_de,gal6_dist_err,gal6_dm_err = distfunc(gal6_p,gal6_ap,grad2,inter2,extin[5],grad_err2,inter_err2)
gal7_dist,gal7_de,gal7_dist_err,gal7_dm_err = distfunc(gal7_p,gal7_ap,grad2,inter2,extin[6],grad_err2,inter_err2)
gal8_dist,gal8_de,gal8_dist_err,gal8_dm_err = distfunc(gal8_p,gal8_ap,grad2,inter2,extin[7],grad_err2,inter_err2)

print(gal_name[0],'=',gal1_dist,'+/-',gal1_de)
print(gal_name[1],'=',gal2_dist,'+/-',gal2_de)
print(gal_name[2],'=',gal3_dist,'+/-',gal3_de)
print(gal_name[3],'=',gal4_dist,'+/-',gal4_de)  # printing out what each galaxy distance is with the error 
print(gal_name[4],'=',gal5_dist,'+/-',gal5_de)
print(gal_name[5],'=',gal6_dist,'+/-',gal6_de)
print(gal_name[6],'=',gal7_dist,'+/-',gal7_de)
print(gal_name[7],'=',gal8_dist,'+/-',gal8_de)

errorfit((2.5e6,6.5e6),gal1_dist_err,11, plt.xlabel('distance error'))  # using the function errorfit to show how the errors on the distance are distributed
 
errorfit((0.25,0.45),gal1_dm_err,7,plt.xlabel('distance modulus error'))

print()

ecg,q = opt.curve_fit(xdata =(gal1_dist_err*1*10**-6), ydata = gal1_dm_err, f=gfunc )
print('Gradient',ecg,'+/-',gerror(q))
ned_gdm_err = ecg * (gal1_dist_err*1*10**-6)

plt.xlim(2.5,6)
plt.plot((gal1_dist_err*1*10**-6),gal1_dm_err, color='red', marker='o', markersize=6, linestyle="none")
plt.ylabel('Distance modulus error')
plt.xlabel('Distance error')               # plotting a scatter graph of distance error and distance modulus error to graphically check if they have any significant correlation 
plt.title('Error Correlation')
plt.plot((gal1_dist_err*1*10**-6),ned_gdm_err, color='blue', marker="None", linewidth=2, linestyle="-")
plt.show()

        
gdist_arr = np.array((gal1_dist,gal2_dist,gal3_dist,gal4_dist,gal5_dist,gal6_dist,\
gal7_dist,gal8_dist),dtype = float)             # creating an array of the distance of each galaxy     

lgdist_arr = np.log10((gdist_arr * 1*10**-6))  # logging the array so its ready to plot in the desired form 

gderr_arr = np.array((gal1_de,gal2_de,gal3_de,gal4_de,gal5_de,gal6_de,gal7_de,gal8_de))  # creating an array of the distance error of each galaxy

pderr_arr = (gderr_arr * 1*10**-6)/(np.log(10)*gdist_arr * 1*10**-6) # propogating the distance error of each galaxy so its in the same form as the distance (logged)

v2_rec = np.log10(v_rec)   # logging the recession velocity array so its ready to plot in the desired form 

gparam, gpcov = opt.curve_fit(f=gfunc, xdata = lgdist_arr, ydata = v2_rec  )# fidnding the parameters for the log(distance) against log(recession velocity)

print()

g_err = gerror(gpcov[0])  
print('Gradient =',gparam,'+/-',g_err)
       # finding error on the parameter 

vrec_new = gparam * lgdist_arr        # finding the expected y values using the optimal parameters, this will be used to plot trend line   

plt.plot(lgdist_arr,v2_rec, color='red', marker='o', markersize=6, linestyle="none")
plt.errorbar(lgdist_arr,v2_rec, xerr=pderr_arr, color='black', linestyle = 'none',capsize = 4)         #  plotting the graph of the log of the galaxy distance against the log of the recession velcoity 
plt.plot(lgdist_arr,vrec_new , color='blue', marker="None", linewidth=2, linestyle="-")
plt.xlabel('log(Galaxy Distance)')
plt.ylabel('log(recession velocity)')
plt.title('Expansion Rate of the Universe')
plt.show()

lgdist_new = v2_rec / gparam

print()
print('X^2 =',chi2(lgdist_arr,pderr_arr,lgdist_new))
print('Reduced X^2 =',chi2(lgdist_arr,pderr_arr,lgdist_new)/(len(lgdist_new) -2)) 


print()

h_con = 10**gparam[0]       # 'unlogging' the gradient and its respective error to get it in the correct forn 
h_con_err = h_con * np.log(10) * g_err 
print('Hubbles constant =',h_con,'+/-',h_con_err) 

print()

u_age = (1/(h_con*1*10**-6))/1000      
u_age_error = (h_con_err*1*10**-6)/(np.square(h_con*1*10**-6)) /1000   # using the hubbles constant and its respective error to work out age of the universe and its error         
print('Age Of The Universe (Gigayears) =',u_age,'+/-',u_age_error)













