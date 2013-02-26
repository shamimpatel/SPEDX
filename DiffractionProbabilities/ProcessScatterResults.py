'''
Processes the data from DiffractionProb.


'''

#0.610865238; #35 degrees
#0.785398163; #45 degrees
#0.34906585; #20 degrees
#0.872664626; #50 degrees
MinAngle = 0.717;
MaxAngle = 0.73;


AngleData = [line.split('\t') for line in open('BraggAngle.txt')];
for l in AngleData:
    del l[-1]; #need to get rid of newline character at the end
    
ProbData = [line.split('\t') for line in open('ScatterProb.txt')];
for l in ProbData:
    del l[-1];

RockingCurveData = [line.split('\t') for line in open('RockingCurve.txt')];
for l in RockingCurveData:
    del l[-1];

KeepPlane = [];

for i in range( 1, len(AngleData[0]) ): 
    KeepPlane.append(False);


for i in range( 1, len(AngleData[0]) ): #cycle from 1 to last part of a line
    #loop down each angledata line
    for AngleLine, ProbLine in zip(AngleData[1:],ProbData[1:]):
        if MaxAngle > float(AngleLine[i]) > MinAngle: # and float(ProbLine[i]) > 1E-5:            
            KeepPlane[i-1] = True;

RemainingPlanes = [];


for bShouldKeep,PlaneName in zip(KeepPlane,AngleData[0][1:] ):
    if(bShouldKeep):
        RemainingPlanes.append( PlaneName );
        
print "Remaining Planes:", RemainingPlanes;
print KeepPlane;

''' 
for i in range( len(AngleData[0])-1, 0, -1 ):
    print i, AngleData[0][i], KeepPlane[i-1];
    if( KeepPlane[i-1] == False):
        for line1, line2 in zip(AngleData,ProbData):
            del line1[i];
            del line2[i];
'''        

Output = open("ProcessedResults.txt" ,'w');

#Output.write( str(len(RemainingPlanes)) + "\t" );
Output.write( "Energy\t" );

for PlaneName in RemainingPlanes:
    Output.write( PlaneName + "\t\t\t");

Output.write("Sum\n");


for angleline, probline, rockingcurveline in zip(AngleData[1:],ProbData[1:],RockingCurveData[1:]):
    Output.write(angleline[0] + "\t"); #print energy for this row
    ProbSum = 0.0;
    for prob in probline[1:]:
        #if MaxAngle > float(angle) > MinAngle: Include all reflections in the sum?
        ProbSum += float(prob); #sum up all probabilities across row
            
    for angle,prob,rockingcurve,i in zip(angleline[1:],probline[1:],rockingcurveline[1:],range(len(probline[1:]))):
        if( KeepPlane[i] ): #only print this line if we're keeping this plane            
            if( ProbSum != 0.0 and MaxAngle > float(angle) > MinAngle ):
                #avoid divide by zero error and also make probability zero if we don't care about the angle
                Output.write( angle + "\t" + str(float(prob)/ProbSum) + "\t" + rockingcurve + "\t");
            else:
                Output.write( angle + "\t" + str(0)                   + "\t" + str(0)       + "\t");
            
    Output.write( str(ProbSum) + "\n");
    
Output.close();

