#monte carlo


#ORBITAL DETERMINATION SSP 2018
#Name: Asher Noel

from math import cos, sin, tan, acos, asin, pi, radians, degrees, isnan
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import warnings
warnings.filterwarnings('ignore')
from f import *
import statistics
import random

name = "2002UXmonte.txt"


##################################################################################################
###################################1) LOAD INPUT FILE, remove colons #############################
with open(name) as data:
    txt = data.read()
    txt = txt.replace(':', " ")

f = open(name, "w")
f.write(txt)
f.close()

file = np.loadtxt(name)

random.seed(90001)
                
a_mc = []
e_mc = []
i_mc = []
O_mc = []
w_mc = []
M_mc = []


def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

def monte_carlo (file):

     #For each observation in the file

    n = 100000
        
    uncertainty_ra1 = radians(file[0, 15]/3600)
    uncertainty_dec1 = radians(file[0, 16]/3600)
    uncertainty_ra2 = radians(file[1, 15]/3600)
    uncertainty_dec2 = radians(file[1, 16]/3600)
    uncertainty_ra3 = radians(file[2, 15]/3600)
    uncertainty_dec3 = radians(file[2, 16]/3600)
        
    ra_var1 = np.random.normal(get_ra_dec_variation(file, 0)[0], uncertainty_ra1, n)
    ra_var2 = np.random.normal(get_ra_dec_variation(file, 1)[0], uncertainty_ra2, n)
    ra_var3 = np.random.normal(get_ra_dec_variation(file, 2)[0], uncertainty_ra3, n)
    dec_var1 = np.random.normal(get_ra_dec_variation(file, 0)[1], uncertainty_dec1, n)
    dec_var2 = np.random.normal(get_ra_dec_variation(file, 1)[1], uncertainty_dec2, n)
    dec_var3 = np.random.normal(get_ra_dec_variation(file, 2)[1], uncertainty_dec3, n)
        
    for sample in range(n):
        method_of_gauss(file,'taylor', ra_var1[sample], dec_var1[sample], ra_var2[sample], dec_var2[sample], ra_var3[sample], dec_var3[sample])
        print(sample)

    a = statistics.median(a_mc)
    e = statistics.median(e_mc)
    i = statistics.median(i_mc)
    O = statistics.median(O_mc)
    w = statistics.median(w_mc)
    M = statistics.median(M_mc)

    print(a, e, i , O, w, M)

    a_array = np.array(a_mc)
    e_array = np.array(e_mc)
    i_array = np.array(i_mc)
    O_array = np.array(O_mc)
    w_array = np.array(w_mc)
    M_array = np.array(M_mc)


    np.save("MCaList", a_array)
    np.save("MCeList", e_array)
    np.save("MCiList", i_array)
    np.save("MCOList", O_array)
    np.save("MCwList", w_array)
    np.save("MCMList", M_array)

    a = np.mean(a_array)
    a_f = reject_outliers(a_array, m=2)
    
    e = np.mean(e_array)
    e_f = reject_outliers(e_array, m=2)


    i = np.mean(i_array)
    i_f = reject_outliers(i_array, m=2)


    O = np.mean(O_array)
    O_f = reject_outliers(O_array, m=2) 


    w = np.mean(w_array)
    w_f = reject_outliers(w_array, m=2)


    M = np.mean(M_array)
    M_f = reject_outliers(M_array, m=2)


    print(a, e, i , O, w, M)

    fig = plt.subplots()
    plt.hist(a_f, alpha = 1, bins = 30 )
    plt.show()

    fig = plt.subplots()
    plt.hist(e_f, alpha = 1, bins = 30)
    plt.show()

    fig = plt.subplots()
    plt.hist(i_f, alpha = 1, bins = 30)
    plt.show()

    fig = plt.subplots()
    plt.hist(O_f, alpha = 1, bins = 30)
    plt.show()

    fig = plt.subplots()
    plt.hist(w_f, alpha = 1, bins = 30)
    plt.show()

    fig = plt.subplots()
    plt.hist(M_f, alpha = 1, bins = 30)
    plt.show()
    

###################################################################################################
####################################2) MAIN FUNCTION ##############################################
def method_of_gauss (file, condition, ra_variation1, dec_variation1, ra_variation2, dec_variation2,ra_variation3, dec_variation3):

    
    ###############################################################################################
    ################################3) Sort Through Global/All Data ###############################
    c = np.array(([0,1,2]))
    JDs = get_JDs (file)
    GTs_all = np.zeros(3)
    k =  0.01720209895
    t_3 = k*(JDs[2] - JDs[1])
    t_1 = k*(JDs[0] - JDs[1])
    GTs_all[0] = t_3 - t_1
    GTs_all[1] = t_1
    GTs_all[2] = t_3
    


    
    RAs, Decs = get_ra_dec_rad(file)

    
    RAs[0] = ra_variation1
    Decs[0] = dec_variation1
    RAs[1] = ra_variation2
    Decs[1] = dec_variation2
    RAs[2] = ra_variation3
    Decs[2] = dec_variation3

    
    rho_hats = get_rho_hats (RAs, Decs)
    Suns = file [:, 12:15]


        #the int() function is used because lists require integer indices, not floats
        
    try: 
        r_2_mag = solve_scaler_eq_lagrange(rho_hats, GTs_all, Suns)[0]
    except IndexError:
        print("No r2 mags!")

                ###################################################################################
                #####################6)First iteration with truncated taylor ######################
    f1, f3, g1, g3 = get_truncated_fg (GTs_all, r_2_mag)
    r1,r2,r3 = get_position_vec(f1, f3, g1, g3, rho_hats, Suns)
    r2_dot = get_velocity_vec(f1,f3,g1,g3,r1,r3)
 

                #The main loop is a for loop because, for 2002 UX, the while loop will never
                #converge, even after 1000+ iterations.


                ###################################################################################
                #####################7)MAIN ITERATION LOOP ########################################
    for iteration in range (50):

                        
        gt,t = correct_light(JDs, f1,f3,g1,g3, rho_hats, Suns, c[1]) 

                    #FOURTH DEGREE TAYLOR SERIES
        if condition == "taylor":
            f1,f3,g1,g3 = taylor_fg (r2, r2_dot, gt)

                    #Update vectors      
        r1,r2,r3 = get_position_vec (f1,f3,g1,g3, rho_hats, Suns)
        r2_dot = get_velocity_vec(f1,f3,g1,g3, r1, r3)

    

            #BREAK CONDITION 1: Test for unreasonably large r2s, sometimes common for 2002 UX. 
    if abs(np.linalg.norm(r2)) > 3 or np.linalg.norm(r2) < 0.1:
        print("    Result: r2 >3")
    
                ###################################################################################
                #####################8)NEW ORBITAL ELEMENTS, EPHEMERIS ############################
                
    elements = OD(r2, r2_dot, JDs[1])
    a_mc.append(elements[0])
    e_mc.append(elements[1])
    i_mc.append(degrees(elements[2]))
    O_mc.append(degrees(elements[3]))
    w_mc.append(degrees(elements[4]))
    M_mc.append(degrees(elements[5]))
                
                

    





#CALLING THE MAIN FUNCTION.
#For the second parameter, "taylor" uses taylor series and "closed" uses closed form functions.         
#method_of_gauss(file, "taylor")
monte_carlo(file)

