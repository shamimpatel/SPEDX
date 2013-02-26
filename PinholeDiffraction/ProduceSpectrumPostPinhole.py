import numpy.random;
import math;

DiffractData = [line.split('\t') for line in open('DiffractResultsPostPinhole.txt')];

print len(DiffractData)

numpy.random.seed(42)

'''

delta(E) = 2.355*w*sqrt( r^2 + F*E/w  )

w = eh pair energy, 3.68ev
r = readout noise (approx 3?)
F = fano factor 0.117
E = xray energy in ev

remove 2.355 for sigma instead of FWHM

'''

Outfile = open('CCDDiffractSpectrum.txt', 'w')

w = 3.65;
r = 3;
F = 0.117;

i = 0;

for XRay in DiffractData:
    Energy = float(XRay[2])*1000; #energy in ev
    if(Energy < 3000.0):
        continue; #throw away secondary fluorescence
    Stdev = w*math.sqrt( r*r + (F*Energy)/w );
    RandomE = numpy.random.normal( Energy, Stdev);
    Outfile.write( XRay[0] + '\t' + XRay[1] + '\t' + str(RandomE/1000.0) + '\n');
    i = i+1;

print i    

Outfile.close();
    
