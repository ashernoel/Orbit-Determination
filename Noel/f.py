#FUNCTIONS FILE
#Name: Asher Noel

from math import cos, sin, tan, acos, asin, pi, radians, degrees, isnan
import numpy as np
import matplotlib.pyplot as plt

name = "2002UXmonte.txt"
with open(name,"r+") as data:
    txt = data.read()
    txt = txt.replace(':', " ")

f = open(name, "w")
f.write(txt)
f.close()

file = np.loadtxt(name)

def convert_ra(ra_input): #Useful for outputting RA Astrometry
    ra = degrees(ra_input)
    ra_h = ra/15
    HH_ra = math.floor(ra_h)
    MM_ra = math.floor(60*(ra_h - HH_ra))
    SS_ra = 60*(60*(ra_h - HH_ra)-MM_ra)
    return str(HH_ra)+ ":" + str(MM_ra) + ":" + str(round(SS_ra,1))

def convert_dec(dec_input): #Useful for outputting DEC Astrometry
    dec = degrees (dec_input)
    if dec > 0: 
        DD_dec = math.floor(dec)
        SS_dec = math.floor(60*(dec - DD_dec))
        AA_dec = 60*(60*(dec - DD_dec)-SS_dec)
    if dec < 0:
        DD_dec = math.floor(-dec)
        SS_dec = math.floor(60*(-dec - DD_dec))
        AA_dec = 60*(60*(-dec - DD_dec)-SS_dec)
        DD_dec = -DD_dec
    return str(DD_dec) + ":" + str(SS_dec) + ":" + str(round(AA_dec,1))


def get_JD (Y, M, D, hh, mm, ss): #Get individual Julian Day
        J_o = 367 * Y - int(7*(Y + int((M + 9)/12))/ 4 ) + int (275* M/9) + D + 1721013.5
        JD = J_o + (hh+mm/60+ss/3600)/24
        return JD

global t_0
t_0 = get_JD(2018,7,22,6,0,0) #EPOCH

def get_JDs(file): #Get Julian Days for the file
    JDs = np.zeros ((np.size(file,0),1))
    for i in range (np.size(file,0)):
        JDs [i,0] = get_JD (file[i,0], file[i,1],file[i,2], file[i,3], file[i,4], file[i,5])
    return JDs

def get_gaussian_time (JDs, c): #convert time to gaussian time 
    GTs = np.zeros((np.size(c,0),3))
    k = 0.01720209895
    for i in range(np.size(c,0)):
        t_3 = k* (JDs[int(c[i,2])] - JDs[int(c[i,1])]) #t3-t2
        t_1 = k* (JDs[int(c[i,0])] - JDs[int(c[i,1])]) #t1-t2
        t = t_3 - t_1
        GTs[i,0] = t
        GTs[i,1] = t_1
        GTs[i,2] = t_3

    return GTs

def get_ra_dec_rad (file): #convert ra, dec to radians  
    RAs = np.zeros ((np.size(file,0)))
    Decs = np.zeros((np.size(file,0)))
    for i in range (np.size(file,0)):
        RAs[i] = radians((file [i, 6] + file [i, 7]/60 + file[i,8]/3600)*15)
        if file[i,9] < 0:
            Decs[i] = radians(file [i, 9] - file [i, 10]/60 - file[i,11]/3600)
        else:
            Decs[i] = radians(file [i, 9] + file [i, 10]/60 + file[i,11]/3600)
    return RAs, Decs

def get_ra_dec_variation (file, i): #For monte carlo, get individual line of Ra, Dec
    RAs = np.zeros (1)
    Decs = np.zeros(1)
    RAs[0] = radians((file [i, 6] + file [i, 7]/60 + file[i,8]/3600)*15)
    if file[i,9] < 0:
        Decs[0] = radians(file [i, 9] - file [i, 10]/60 - file[i,11]/3600)
    else:
        Decs[0] = radians(file [i, 9] + file [i, 10]/60 + file[i,11]/3600)
    return RAs, Decs

def get_rho_hats (ra, dec): 
    rhos = np.zeros ((np.size(file, 0),3))
    for i in range (np.size(file,0)):
        rhos [i,0] = cos (ra[i])*cos(dec[i])
        rhos [i,1] = sin (ra[i])*cos(dec[i])
        rhos [i,2] = sin (dec[i])
    return rhos

def solve_scaler_eq_lagrange(rho, GT, sun):

    mu = 1

    A_1 = GT[2]/GT[0]
    B_1 = A_1 / 6 * (GT[0]**2 - GT[2]**2)
    A_3 = - GT[1]/GT[0]
    B_3 = A_3/6 * (GT[0]**2 - GT[1]**2)

    
    D_0 = np.dot(rho[0], np.cross (rho[1], rho[2]))
    D_11 = np.dot(np.cross(sun[0],rho[1]), rho[2])
    D_12 = np.dot(np.cross(sun[1],rho[1]), rho[2])
    D_13 = np.dot(np.cross(sun[2],rho[1]), rho[2])
    D_21 = np.dot(np.cross(rho[0],sun[0]), rho[2])
    D_22 = np.dot(np.cross(rho[0],sun[1]), rho[2])
    D_23 = np.dot(np.cross(rho[0],sun[2]), rho[2])
    D_31 = np.dot(rho[0],np.cross(rho[1], sun[0]))
    D_32 = np.dot(rho[0],np.cross(rho[1], sun[1]))
    D_33 = np.dot(rho[0],np.cross(rho[1], sun[2]))
    

    A = (A_1 * D_21 - D_22 + A_3*D_23) / (-D_0)
    B = (B_1 * D_21 + B_3 * D_23)/ (-D_0)

    E = -2 *(np.dot (rho[1], sun[1]))
    F = (np.linalg.norm(sun[1]))**2
    
    a = - (A**2 + A*E+F)
    b = - mu *( 2*A*B + B*E)
    c = - mu**2 * B **2 

    fig = plt.subplots()
    x = np.linspace(0.2, 1.9, 100)
    plt.plot(x, x**8 + a*x**6 + b**3 + c, label = "Eqn of Lagrange")
    plt.xlabel("2002 UX distance [AU]", fontsize = 14, fontname = "Times New Roman")
    plt.ylabel("f(r)", fontsize = 14, fontname = "Timew New Roman")
    plt.title("Scalar Equation of Lagrange", fontsize = 14, fontname = "Times New Roman")
    plt.show()
    
    r = np.roots([1,0,a,0,0,b,0,0,c])
    r = r[np.isreal(r)].real
    r = r[r > 0]

    #Get rid of -rhos and complex roots
    for i in range(len(r)):
        rho = A + B/r[i]**(3)
        if rho < 0:
            r[i] = -69
        if r[i] > 10:
            r[i] = -69
    r = r[r>0]
    return r

def get_truncated_fg (GT, r2):
    mu = 1
    f_1 = 1 - (mu*GT[1]**2)/(2*r2**3)
    f_3 = 1 - (mu*GT[2]**2)/(2*r2**3)
    g_1 = GT[1] - (mu*GT[1]**3)/(6*r2**3)
    g_3 = GT[2] - (mu*GT[2]**3)/(6*r2**3)
    return f_1, f_3, g_1, g_3

def get_position_vec(f1,f3,g1,g3,rho,sun):

    c1 = (g3)/(f1*g3-g1*f3)
    c3 = (-g1)/(f1*g3-g1*f3)
    D_0 = np.dot(rho[0], np.cross (rho[1], rho[2]))
    D_11 = np.dot(np.cross(sun[0],rho[1]), rho[2])
    D_12 = np.dot(np.cross(sun[1],rho[1]), rho[2])
    D_13 = np.dot(np.cross(sun[2],rho[1]), rho[2])
    D_21 = np.dot(np.cross(rho[0],sun[0]), rho[2])
    D_22 = np.dot(np.cross(rho[0],sun[1]), rho[2])
    D_23 = np.dot(np.cross(rho[0],sun[2]), rho[2])
    D_31 = np.dot(rho[0],np.cross(rho[1], sun[0]))
    D_32 = np.dot(rho[0],np.cross(rho[1], sun[1]))
    D_33 = np.dot(rho[0],np.cross(rho[1], sun[2]))

    c2 = -1

    rho_1 = (c1*D_11 + c2*D_12 + c3*D_13)/(c1*D_0)
    rho_2 = (c1*D_21 + c2*D_22 + c3*D_23)/(c2*D_0)
    rho_3 = (c1*D_31 + c2*D_32 + c3*D_33)/(c3*D_0)
    r1 = np.zeros(3)
    r2 = np.zeros(3)
    r3 = np.zeros(3)

    for i in range(3):

        r1 [[i]] = rho_1 * rho[0,i] - sun[0,i]
        r2 [[i]] = rho_2 * rho[1,i] - sun[1,i]
        r3 [[i]] = rho_3 * rho[2,i] - sun[2,i]
    return r1,r2,r3

def get_velocity_vec(f1,f3,g1,g3, r1,r3):
    d_1 = -(f3)/(f1*g3-f3*g1)
    d_3 = f1/(f1*g3-f3*g1)
    
    r2_dot = np.zeros(3)
    for i in range(3):
        r2_dot[i]=d_1*r1[i]+d_3*r3[i]
        
    return r2_dot

def correct_light (JDs, f1,f3,g1,g3, rho, sun, c):
    c1 = (g3)/(f1*g3-g1*f3)
    c3 = (-g1)/(f1*g3-g1*f3)
   
    D_0 = np.dot(rho[0], np.cross (rho[1], rho[2]))
    D_11 = np.dot(np.cross(sun[0],rho[1]), rho[2])
    D_12 = np.dot(np.cross(sun[1],rho[1]), rho[2])
    D_13 = np.dot(np.cross(sun[2],rho[1]), rho[2])
    D_21 = np.dot(np.cross(rho[0],sun[0]), rho[2])
    D_22 = np.dot(np.cross(rho[0],sun[1]), rho[2])
    D_23 = np.dot(np.cross(rho[0],sun[2]), rho[2])
    D_31 = np.dot(rho[0],np.cross(rho[1], sun[0]))
    D_32 = np.dot(rho[0],np.cross(rho[1], sun[1]))
    D_33 = np.dot(rho[0],np.cross(rho[1], sun[2]))

    c2 = -1

    rho_1 = (c1*D_11 + c2*D_12 + c3*D_13)/(c1*D_0)
    rho_2 = (c1*D_21 + c2*D_22 + c3*D_23)/(c2*D_0)
    rho_3 = (c1*D_31 + c2*D_32 + c3*D_33)/(c3*D_0)
    
    #c = 173.145 AU/day for speed of light
    t = np.zeros(3)
    t[[0]] = JDs[0] - rho_1/173.145
    t[[1]] = JDs[1] - rho_2/173.145
    t[[2]] = JDs[2] - rho_3/173.145

    gts = np.zeros(3)
    k = 0.01720209895

    gts[2] = k* (t[2] - t[1])
    gts[1] = k* (t[0] - t[1])
    gts[0] = gts[2]-gts[1]
    return gts ,t #Gaussian time and julian day output

        
def get_dE (gt, a, r, r_dot,e): #Broken function for closed form 

    n = (1/a**(3))**(1/2)
    d = (r * r_dot )/ (n*a**2)

    x_n =n*( gt )
    x = 1
    
    while abs(x_n - x) > 1e-12:
            
            x = x_n
            f_x = x - (1 - r/a) * sin(x)+ d *(1-cos(x)) - n*(gt)
            f1_x = 1 - (1 - r/a) * cos(x) + d * sin(x)
            x_n = x - f_x/f1_x
            
    return x_n


def closed_fg(r, r_dot, gt, t): #Another broken function, don't try it

    a = (2/np.linalg.norm(r) - (np.linalg.norm(r_dot)**2))**(-1)
    e = (1 - (np.linalg.norm(np.cross(r.transpose(), r_dot.transpose())))**(2)/a)**(1/2)
    
    n = (1/a**(3))**(1/2)
    print (n*gt[2])
    
    E1 = get_dE(gt[1], a, np.linalg.norm(r), np.linalg.norm(r_dot),e)
    E3 = get_dE(gt[2], a, np.linalg.norm(r), np.linalg.norm(r_dot),e)
    
    print(E1, E3)

    f1 = 1 - (a /np.linalg.norm(r)) * (1-cos (E1))
    f3 = 1 - (a /np.linalg.norm(r)) * (1-cos (E3))
    g1 = gt[1] + 1/n * (sin (E1)-E1)
    g3 = gt[2] + 1/n * (sin (E3)-E3)

    return f1,f3,g1,g3

def taylor_fg (r2, r2_dot, gt): 
    u = 1/np.linalg.norm(r2)**3
    z = np.dot (r2, r2_dot)/np.linalg.norm(r2)**2
    q = np.dot(r2_dot, r2_dot)/np.linalg.norm(r2)**2 - u
   
    
    f1 = 1 - gt[1]**2*u/2 + 1/2*u*z*gt[1]**3 + 1/24*(3*u*q-15*u*z**2 + u**2)*gt[1]**4
    f3 = 1 - gt[2]**2*u/2 + 1/2*u*z*gt[2]**3 + 1/24*(3*u*q-15*u*z**2 + u**2)*gt[2]**4
    
    g1 = gt[1] - 1/6*u*gt[1]**3 + 1/4*u*z*gt[1]**4
    g3 = gt[2] - 1/6*u*gt[2]**3 + 1/4*u*z*gt[2]**4
    return f1, f3, g1, g3

def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)
    if cosine < 0 and sine > 0: #2
        return acos(cosine)
    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)
    if cosine > 0 and sine < 0: #4
        return 2*pi +asin(sine)
    
def check_radians(rad):
    while rad < 0:
        rad += 2*pi
    while rad > 2*pi:
        rad -= 2*pi
    return rad

def OD(pos_eq, vel_eq, JD):

    #Rotate to ecliptic
    pos_ec1 = eq_to_ec(pos_eq)[[0]]
    vel_ec1 = eq_to_ec(vel_eq)[[0]]

    pos_ec = np.zeros(3)
    vel_ec = np.zeros(3)
    for i in range (3):
        pos_ec [i] = pos_ec1[0,i]
        vel_ec [i] = vel_ec1[0,i]
        
    #1: Semimajor Axis (a)
    r = np.linalg.norm(pos_ec)
    v = np.linalg.norm(vel_ec)
    a = (2/r - v**2)**(-1)

    #2: Eccentricity (e)
    e = (1- (np.linalg.norm(np.cross(pos_ec.transpose(), vel_ec.transpose())))**(2)/a)**(1/2)

    #3: Inclination (i)
    k = 0.01720209894
    h = k * (np.cross(pos_ec.transpose(), vel_ec.transpose()))
    i = acos(h[2]/np.linalg.norm(h)) #radians inclination

    #4: Omega (O)
    sinO = h[0]/(np.linalg.norm(h)*sin(i))
    cosO = -h[1]/(np.linalg.norm(h)*sin(i))
    Omega = findQuadrant (sinO, cosO)#radians longitude of the ascending node Omega

    
    #5: omega, mean anomaly
    sinf = (np.dot(pos_ec.transpose(), vel_ec))/(e*r) * (a*(1-e**2))**(1/2)
    cosf = 1/e * (a*(1-e**2)/r -1)
    f = findQuadrant (sinf, cosf)

    sinwf = pos_ec[2]/(r*sin(i))
    coswf = (cos(Omega))**(-1)* (pos_ec[0]/r + cos(i) * sinwf * sin(Omega))
    w = check_radians(findQuadrant(sinwf, coswf)-f) #Radians mean anomaly

    n = k*(1/a**(3))**(1/2)
    t_2 = JD

    #f, E in same plane, so Lawson said, "It's a quadrant rule" 
    if f > pi:
        E = -acos (1/e* (1- r/a))
    else:
        E = acos (1/e* (1- r/a))
        
    M_2 = E - e * sin (E)
    M_0 = check_radians(M_2 - n *(t_2 - t_0))    

    return a, e, i, Omega, w, M_0

def eq_to_ec(eq): #Rotation
    e = radians(23.4352)
    ec = np.zeros(3)
    eq_trans = np.matrix(([1,0,0],[0,cos(e), -sin(e)],[0,sin(e),cos(e)]))
    eq_inv = np.linalg.inv(eq_trans)
    ec = np.dot(eq_inv, eq)

    return ec

def ec_to_eq(ec): #Rotation
    e = radians(23.4352)
    eq = np.zeros((3,1))
    eq_trans = np.matrix(([1,0,0],[0,cos(e), -sin(e)],[0,sin(e),cos(e)]))
    eq = np.dot(eq_trans, ec)
    return eq

def get_M (elements , t,i): #Get M and add to elements array
    k = 0.01720209895

    n = k*(1/elements[1]**(3))**(1/2)
    
    M = n*(t[i] - t_0) + elements[5]
    elements[6] =  M
    return None
        
def get_E2 (elements): #Get E using Kepler Equationa and add to elements array
    e = elements[0]
    M = elements[6]
    E_g = M
    M_g = E_g - e*sin(E_g)
    M_true = M
    while abs(M_g - M_true) > 10**(-10):
        M_g = E_g - e*sin(E_g)
        E_g = E_g - (M_true - (E_g - e*sin(E_g))) / ( e*cos(E_g)-1)
    elements[7] = E_g
    return None

def get_cartesian(orb): #Rotation
    xyz = np.zeros((3,1))
    xyz[0, 0] = orb[1]*cos(orb[7])-orb[1]*orb[0]
    xyz[1, 0] = orb[1]*(1-orb[0]**2)**(1/2)*sin(orb[7])
    xyz[2, 0] = 0
    return xyz

def ecliptic(orb, cartesian): #Rotation to ecliptic
    w = orb[4] 
    i = orb[2]
    O = orb[3]
    r1 = np.matrix(([cos(w), -sin(w), 0], [sin(w), cos(w), 0], [0,0,1]))
    r2 = np.matrix(([1,0,0],[0,cos(i),-sin(i)],[0,sin(i),cos(i)]))
    r3 = np.matrix(([cos(O), -sin(O), 0],[sin(O), cos(O), 0], [0,0,1]))
    ec = np.zeros((3,1))
    ec = np.dot(r3, np.dot(r2, np.dot(r1, cartesian)))
    return ec


def ephemeris(elements, sun,t,k):

    #Put all elements into one array
    orbital_elements = np.array([elements[1], elements[0], elements[2], elements[3], elements[4], elements[5],0,0])

    get_M(orbital_elements, t,k) #k will be iterated to get local time with t in JD
    get_E2(orbital_elements)

    ct_coord = get_cartesian(orbital_elements)
    ec_coord = ecliptic(orbital_elements, ct_coord)
    eq_coord = ec_to_eq(ec_coord)
    
    es_eq_coord = np.zeros((3,1))
    for j in range (3):
        es_eq_coord [j,0] = sun[k,j] #iteration specific sun vector
    
    rho = eq_coord + es_eq_coord 
    rho_hat = rho / np.linalg.norm(rho)
    dec = asin(rho_hat[2])
    a1 = (rho_hat[0]/cos(dec))
    a2 = (rho_hat[1]/cos(dec))
    asc = findQuadrant(a2, a1)
    return asc, dec

def coord (t, j, vecs, delta, sun): #EPEHEMERIS generator for the partial derivative, returning ras and decs

    #convert 6,1 "vecs" matrix back to r2, r2_dot with slight iteration dependent "j" adjustment by delta
    r2 = np.zeros(3)
    r2_dot = np.zeros(3)
    
    for i in range (3):
        vecs [j] += delta 
        r2[i] = vecs [i]
        r2_dot [i] = vecs [i+3]

    #Get RAs and Decs for the first half of the partial derivative.     
    elements = OD(r2, r2_dot, t[1])
    ra1, dec1 = np.zeros(5), np.zeros(5)
    for k in range (5):
        ra1[k], dec1[k] = ephemeris(elements, sun, t, k)


    #Convert to new r_2, r2_dot
    for i in range (3):
        vecs [j] -= 2 * delta  #2*delta to make up for the first + delta 
        r2[i] = vecs [i]
        r2_dot [i] = vecs [i+3]
        
    elements = OD(r2, r2_dot, t[1])
    ra2, dec2 = np.zeros(5), np.zeros(5)
    for k in range (5):
        ra2[k], dec2[k] = ephemeris(elements, sun, t,k)
        
    return ra1, dec1, ra2, dec2


def partial (t, vecs, i, sun): #The partial derivative function
    delta = 1e-3
    
    p_ra = 0
    p_dec = 0
    
    p_ra = np.sum(coord(t, i, vecs, delta, sun)[0] - coord(t, i, vecs, delta, sun)[2])/(2*delta)
    p_dec = np.sum(coord(t, i, vecs, delta, sun)[1] - coord(t, i, vecs, delta, sun)[3])/(2*delta)
    
    return p_ra, p_dec

def differential_correction (r, r_dot, ra, dec, RAs, Decs, t, sun):

    #Combine r and r_dot into one 6 item vector

    vecs = np.zeros (6)
    for i in range (3):
        vecs [i] = r [i]
        vecs [i+3] = r_dot [i]

    resids = np.zeros((6,1))
    p_trans = np.zeros((6,6))

    oc_ra = np.zeros(5)
    oc_dec = np.zeros(5)
    for i in range(5):
        oc_ra[i] = RAs[i] - ra[i]
        oc_dec[i] = Decs[i] - dec[i]

    res = np.sum(oc_ra + oc_dec) #The scalar sum of the O-C values. 

    p=np.zeros(6) #Will act as 6,1 matrix of six summed partial derivatives

    for i in range (6):
        p[i] += partial(t, vecs, i, sun)[0] #returns RA partial
        p[i] += partial(t, vecs, i, sun)[1] #returns DEC partial
        
    for i in range (6):
        for j in range(6):
            p_trans [i,j] = p[i] * p[j] #Creates 6x6 J transformation matrix
        resids [i,0] = res * p[i] #Creates 6,1 a matrix

    correction = np.zeros((6,1))
    correction = np.dot(np.linalg.pinv(p_trans),resids) 
    #inv = np.linalg.inv(p_trans) #This results in a singular matrix
    #correction = np.linalg.lstsq(p_trans, resids)[0]

    return correction

def get_RMS (elements_0, elements, ra0, dec0, Suns, JD_correct):

    #Check for number of observations inputted through the sun vector parameter
    if np.size(Suns,0) == 5:
        obs = 5
    else:
        obs = 3


    ra1, dec1 = np.zeros(obs), np.zeros(obs)
    ra2, dec2 = np.zeros(obs), np.zeros(obs)
    for k in range (obs):
        ra1[k], dec1[k] = ephemeris(elements_0, Suns, JD_correct, k)
        ra2[k], dec2[k] = ephemeris(elements, Suns, JD_correct, k)

    res1 = 0
    res2 = 0
    for i in range (np.size(ra1)):
        res1 += (ra0[i] - ra1[i])**2
        res1 += (dec0[i] - dec1[i])**2
        res2 += (ra0[i] - ra2[i])**2
        res2 += (dec0[i] - dec2[i])**2
        
    RMS1 = (res1/(5*2-6))**(1/2)
    RMS2 = (res2/(5*2-6))**(1/2)

    return RMS1, RMS2
