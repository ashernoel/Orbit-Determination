
#ORBITAL DETERMINATION SSP 2018
#Name: Asher Noel

from math import cos, sin, tan, acos, asin, pi, radians, degrees, isnan
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from f2 import *


###################################################
##################INPUT FILE HERE##################
name = "2002UX.txt"



##################################################################################################
###################################1) LOAD INPUT FILE, remove colons #############################
with open(name) as data:
    txt = data.read()
    txt = txt.replace(':', " ")

f = open(name, "w")
f.write(txt)
f.close()

file = np.loadtxt(name)


###################################################################################################
####################################2) MAIN FUNCTION ##############################################
def method_of_gauss (file, condition):

    
    ###############################################################################################
    ################################3) Sort Through Global/All Data ###############################

    if np.size(file,0) == 5:
        c = np.array(([0,1,4], [0,2,4],[0,3,4]))
    elif np.size(file,0) == 4:
        c = np.array(([0,1,3], [0,2,3]))
    else:
        c = np.array(([1,2,3]))
    JDs = get_JDs (file)
    GTs_all = get_gaussian_time (JDs, c)
    RAs, Decs = get_ra_dec_rad(file)
    rho_hats = get_rho_hats (RAs, Decs)
    Suns = file [:, 12:15]

    #Create empty matrices that will be filed during the iterations and used for Diff Correction
    r2_all = np.zeros ((np.size(c,0), 3))
    r2_dot_all = np.zeros ((np.size(c,0), 3))
    ra_c = np.zeros((np.size(c,0), np.size(file,0)))
    dec_c = np.zeros((np.size(c,0), np.size(file,0)))
    JD_correct = np.zeros((np.size(file,0),1))


    ###############################################################################################
    ################################4) Loop Through Combinations ##################################
    for i in range (np.size(c,0)):


        ###########################################################################################
        ############################5) Create Combination Specific Data Structures ################
        print("\nCombination: " + str(i+1) + " " + str(c[i]))
        local_GTs = GTs_all[i]
        local_RAs = RAs[[[int(c[i,0]), int(c[i,1]), int(c[i,2])]]]
        local_Decs = Decs[[[int(c[i,0]), int(c[i,1]), int(c[i,2])]]]
        local_rho_hats = rho_hats[[[int(c[i,0]), int(c[i,1]), int(c[i,2])]]]
        local_suns = file[[int(c[i,0]), int(c[i,1]), int(c[i,2])],12:15]
        local_times = JDs[[[int(c[i,0]), int(c[i,1]), int(c[i,2])]]]

        #the int() function is used because lists require integer indices, not floats
        
        r_2_mag = solve_scaler_eq_lagrange(local_rho_hats, local_GTs, local_suns)
        
        #Loop through every positive rho, real root (of which there is usually 1)
        for j in range (np.size(r_2_mag)):


                ###################################################################################
                #####################6)First iteration with truncated taylor ######################
                print("  r2: ", r_2_mag[j])
                f1, f3, g1, g3 = get_truncated_fg (local_GTs, r_2_mag[j])
                r1,r2,r3 = get_position_vec(f1, f3, g1, g3, local_rho_hats, local_suns)
                r2_dot = get_velocity_vec(f1,f3,g1,g3,r1,r3)
                print("    1st Iteration: \n        r2_vec: " + str(r2) + " AU \n        r2_dot: "
                      + str(0.01720209895*r2_dot)+ " AU/day")

                #The main loop is a for loop because, for 2002 UX, the while loop will never
                #converge, even after 1000+ iterations.


                ###################################################################################
                #####################7)MAIN ITERATION LOOP ########################################
                for iteration in range (100):

                        
                    gt,t = correct_light(local_times, f1,f3,g1,g3, local_rho_hats, local_suns, c[i]) 

                    for k in range (3):
                        JD_correct [c[i,k],0] = t[k] #Update JDs for Differential Correction
                    #CLOSED FORM FUNCTIONS ARE BROKEN
                    if condition == "closed":
                        f1,f3,g1,g3 = closed_fg (r2, r2_dot, gt, t)

                    #FOURTH DEGREE TAYLOR SERIES
                    if condition == "taylor":
                        f1,f3,g1,g3 = taylor_fg (r2, r2_dot, gt)

                    #Update vectors      
                    r1,r2,r3 = get_position_vec (f1,f3,g1,g3, local_rho_hats, local_suns)
                    r2_dot = get_velocity_vec(f1,f3,g1,g3, r1, r3)

                #BREAK CONDITION 1: Test for NaN, a common output of the closed form
                if isnan(r2[0]):
                    print("    Result: NAN")
                    break

                #BREAK CONDITION 2: Test for unreasonably large r2s, sometimes common for 2002 UX. 
                if np.linalg.norm(r2) > 5:
                    print("    Result: r2 >5")
                    break

                print("    100th Iteration: \n        r2_vec: " + str(r2) +
                      " AU\n        r2_dot: " + str(0.01720209895*r2_dot) + " AU/day")
                print("         Range: " + str(np.linalg.norm(r2)) + " AU")


                ###################################################################################
                #####################8)NEW ORBITAL ELEMENTS, EPHEMERIS ############################
                
                elements = OD(r2, r2_dot, JDs[int(c[i,1])])
                
                ra, dec = np.zeros(3), np.zeros(3)
                for k in range (3):
                    #Define local RA, DEC
                    ra[k], dec[k] = ephemeris(elements, local_suns, t,k)

                    #Update global RA, DEC
                    ra_c [i, c[i,k]] = ra[k]
                    dec_c [i, c[i,k]] = dec[k]

                    #Update global r2, r2_dot
                    r2_all [i, k] = r2 [k]
                    r2_dot_all [i, k] = r2_dot [k]

        

                #Print outputs
                print("            a = " + str(elements[0]) + "AU \n            e = " +
                      str(elements[1]) + "\n            i = " + str(degrees(elements[2])) +
                      " deg \n            O = " + str(degrees(elements[3])) +
                      " deg \n            w = " + str(degrees(elements[4])) +
                      " deg \n            M = " + str(degrees(elements[5])) +
                      " deg at JD = " + str(t[1]))
                print("    RA_c: " + str(ra) + " Dec_c: " + str(dec))
                print("    RA_o: " + str(local_RAs) + " Dec_o: " + str(local_Decs))

        
    ###############################################################################################
    ################################9) Differential Correction ####################################
                
    #Choose the row by hand that will serve as the "average" vector for Differential Correction.
    #By Default, this is the first r2, r2_dot
    r2_av = r2_all [0,:]
    r2_dot_av = r2_dot_all [0,:]

    #Calculate orbital elements and RA, Decs for all FIVE observations
    elements_old = OD(r2_av, r2_dot_av, JD_correct[1,0])

    print(JD_correct)
    print(JDs)
    print(Suns)
    ra1, dec1 = np.zeros(5), np.zeros(5)
    for k in range (5):
        ra1[k], dec1[k] = ephemeris(elements_old, Suns, JDs, k)
    
    print(ra1, dec1)
    print(RAs, Decs)
    #Return the (6,1) correction (x) matrix, and then adjust r2, r2_dot
    correction = differential_correction (r2_av, r2_dot_av, ra1, dec1, RAs, Decs, JDs, Suns)

    for k in range (3):
        r2_av[k] = +correction [k] + r2_av[k]
        r2_dot_av [k] = +correction [k+3] + r2_dot_av[k]

    #FINAL ORBITAL ELEMENT CALCULATION (with RMS)
    elements = OD(r2_av, r2_dot_av, JD_correct[1,0])
    RMS = get_RMS(elements_old, elements, RAs, Decs, Suns, JDs)

    
    print(RMS)  
    print("            a = " + str(elements[0]) + "AU \n            e = " + str(elements[1]) +
          "\n            i = " + str(degrees(elements[2])) + " deg \n            O = " +
          str(degrees(elements[3])) + " deg \n            w = " + str(degrees(elements[4])) +
          " deg \n            M = " + str(degrees(elements[5])) +
          " deg at JD = " + str(get_JD(2018,7,22,6,0,0)))
                

#CALLING THE MAIN FUNCTION.
#For the second parameter, "taylor" uses taylor series and "closed" uses closed form functions.         
method_of_gauss(file, "taylor")
