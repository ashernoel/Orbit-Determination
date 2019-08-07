from vpython import *
import numpy as np

#Get the classical orbital elements and determine M at ephemeris time:

#Defining Old Asteroiod
a = 1.523653749051593E+00
e = 9.332582302262121E-02
M = radians(3.303898234525319E+02)
Oprime = radians(4.950665272307243E+01)
iprime = radians( 1.848148620293442E+00)
wprime = radians(2.866558291231743E+02)

#Defining Earth
a_earth = 9.996724564124407E-01
e_earth = 1.710413998386059E-02
M_earth = radians(1.972756360820545E+02)
Oprime_earth = radians(1.862973760375422E+02)
iprime_earth = radians(2.037179991363645E-03)
wprime_earth = radians(2.759112299409064E+02)
p_earth = sqrt (4*pi**(2) *a_earth**(3))
r_earth = vector (0,0,0)
earth = sphere (pos = r_earth * 150, radius = (15), color = color.blue)
earth.trail = curve (color=color.blue)
M_t_earth = 2*pi/p_earth * 0 + M_earth

#Defining 2002 UX
a_UX = 1.4735853
e_UX = 0.1633837
M_UX = radians(253)
Oprime_UX = radians (84)
iprime_UX = radians (20)
wprime_UX = radians (263)
p_UX = sqrt (4*pi**(2) *a_UX**(3))
r_UX = vector (0,0,0)
UX = sphere (pos = r_UX * 150, radius = (7), color = color.red)
UX.trail = curve (color=color.red)
M_t_UX = 2 * pi / p_UX * 0 + M_UX

#Get_E function for all orbits
def get_E (M):

    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = M
    while abs(Mguess - Mtrue) > 1e-004:
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (Mtrue - (Eguess - e*sin(Eguess))) / ( e*cos(Eguess)-1)
    return Eguess

#Transformation matrix for all orbits
def transform(E, a, e, i, w, Omega):
                   
    cartesian = vector(a*cos(E)- a*e, a*sqrt(1-e**2)*sin(E), 0)
    v1 = rotate(cartesian, angle = -w, axis = vector(0, 0, 1))
    v2 = rotate(v1, angle = i, axis = vector(1, 0, 0))
    v3 = rotate(v2, angle = Omega, axis = vector(0, 0, 1))

    return v3

#Leftover definitions of asteroid and beginning of its orbit

time = 0
dt = .05

period = sqrt (4*pi**(2) *a**(3))

r1ecliptic = vector (0,0,0)
asteroid = sphere(pos = r1ecliptic * 150, radius=(10), color=color.white)
asteroid.trail = curve(color=color.white)
M_true = 2*pi/period * (time) + M
E_true = get_E (M_true)
v3 = transform(E_true, a, e, iprime, wprime, Oprime)
asteroid.pos = v3 * 150


sun = sphere(position= vector(0,0,0), radius=(30), color=color.yellow)

#Placing of the earth and Labeling
E_t_earth = get_E (M_t_earth)
v_earth = transform (E_t_earth, a_earth, e_earth, iprime_earth, wprime_earth, Oprime_earth)
earth.pos = v_earth * 150

#Placing of 2002 UX and Labeling
E_t_UX = get_E (M_t_UX)
v_UX = transform (E_t_UX, a_UX, e_UX, iprime_UX, wprime_UX, Oprime_UX)
UX.pos = v_UX * 150


while (1 == 1):

    #Update Asteroid Position
    rate(10)
    time += dt

    
    M_true = 2*pi/period * (time) + M
    E_true = get_E (M_true)
    v3 = transform(E_true, a, e, iprime, wprime, Oprime)
    asteroid.pos = v3 * 150
    asteroid.trail.append(pos = asteroid.pos)
    
    
    #Update Earth Position
    M_t_earth = 2*pi/p_earth * time + M_earth
    E_t_earth = get_E (M_t_earth)
    v_earth = transform (E_t_earth, a_earth, e_earth, iprime_earth, wprime_earth, Oprime_earth)
    earth.pos = v_earth * 150
    earth.trail.append( pos = earth.pos)

    #Update 2002 UX Position
    M_t_UX = 2 * pi / p_UX * time + M_UX
    E_t_UX = get_E (M_t_UX)
    v_UX = transform (E_t_UX, a_UX, e_UX, iprime_UX, wprime_UX, Oprime_UX)
    UX.pos = v_UX * 150
    UX.trail.append(pos = UX.pos)



    
    
    
    
