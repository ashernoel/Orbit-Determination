# Orbital-Determination
This program arrives at differentially corrected orbital elements that describe and can be used to model the motion of an asteroid after image reduction and processing using least-squares plate reduction astrometry. 

This code was developed as a part of the Summer Science Program in Astrophysics at the New Mexico Institute of Technology during the summer of 2018. 

## Output

<img src="/Images/observation1.png" alt="Default Login Screen" width="250"/>  <img src="/Images/observation2.png" alt="Default Login Screen" width="250"/> 
<img src="/Images/montecarlo.png" alt="Default Login Screen" width="200"/>  <img src="/Images/orbitvisualization.png" alt="Default Login Screen" width="300"/> 

| Orbital Element        |   Output           |  JPL    |
|:-------------:|:-------------:|:--------:|
| *a*  | 1.4692 AU | 1.4734 AU|
| *e*     | 0.1626      |   0.1633 |
| *i* | 20.07&deg;       |  20.208°|
| *&Omega;* | 263.83°  |  263.90°|
| *&omega;* | 84.739°   |  84.273°|
| *M* | 319.56° |  319.94°|

## Files:

### NoelOD.py  &  f2.py

This is the main function. It imports functions from f2.py. 

To run it, 
1.  Open NoelOD.py and f2.py;
2.  Replace "input file here.txt" with the input file;
3. Save f2.py and run NoelOD.py; and 
4. Debug.

### MonteCarlo.py  &  f.py

This is a modified Method of Gauss code that works with MonteCarlo and uses basic functions from f.py.
Replace the "input file here.txt" with the input file in both .py files to run. The input file can only have 3 observations. Uncertainties should be in the final two columns. 

### OrbitVisualization.py

2002 UX's orbit visualized with Mars, Earth, and the Sun.
