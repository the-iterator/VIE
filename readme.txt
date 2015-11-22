--------------------------------------------------------------------------------

Volume Integral Equation (VIE) solver 

viepy2.py works in Python 2.7
viepy3.py works in both Python 2.7 and Python 3.*

Neil Budko (c) 2015

--------------------------------------------------------------------------------

What it does:

Given the source (location and polarization of a point-dipole 
or direction and polarization of a plane wave), the geometry, and 
the permittivity of an inhomogeneous scatterer, this Python 2.7/3 
code computes the complex amplitude of the monochromatic electric 
field over the computational domain (three components), displays 
it with slices, and saves it as Matlab arrays

--------------------------------------------------------------------------------

Requires python 2.7/3.* with modules:

scipy
matplotlib
time
sys
os
inspect

also recommended: ipython or ipython3

--------------------------------------------------------------------------------

Main folder (e.g. ~/VIE/) should contain the following files

VIEforce.py
incident.py
addobject.py
green.py
matvec.py
visualization.py

subfolder (VIE/Input) with at least one *.inp input file

subfolder (VIE/Output) may be empty, will contain *.mat output files

--------------------------------------------------------------------------------

To run the code:

1. Create your own or choose an existing input file inside /Input/ folder. 
Set its name inside viepy2.py/viepy3.py

2. Go to the main folder (/VIE/) and start ipython in terminal
(if you are not in the main folder you can "cd" to it from inside ipython)

3. In ipython type: run viepy2.py input_file 
where 'input_file' is the name of your *.inp file without extension 
(in ipython3 type: run viepy3.py input_file)

4. Watch the text output and wait for four figures to appear

5. When done, close figures. You can change the input file and run 
viepy2.py/viepy3.py again, or exit ipython/ipython3 by typing 'exit'.

5. For 3D visualization use Matlab with the data arrays in /Output/ folder.
For example, run the matslice.m script which is inside the /Output/ folder
--------------------------------------------------------------------------------

