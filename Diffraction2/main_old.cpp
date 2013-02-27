#include <iostream>
#include "Vector.h"
#include "math.h"
#include "FormFactorData.h"
#include "LatticePlane.h"
using namespace std;


float UpperLimit( float Value, float Limit)
{
  if( Value > Limit )
  {
    return Limit;
  }
  else
  {
    return Value;
  }
}

bool ApproxEqual( float a, float b, float tolerance)
{
  return fabs( a - b ) <= tolerance;
}


int main ()
{  

  FormFactorData IronFormFactor( 0.0f, 1.0f, 5000);
  IronFormFactor.LoadData("FeFormFactor.txt");
  
  float a0 = 5; //lattice constant 5 angstroms

  Vector a1(1,0,0);
  Vector a2(0,1,0);
  Vector a3(0,0,1);

  a1 = a1*a0;
  a2 = a2*a0;
  a3 = a3*a0;

  float UnitCellVol = a2.Cross(a3).Dot(a1);
  
  Vector b1 = (a2.Cross(a3))/UnitCellVol;
  Vector b2 = (a3.Cross(a1))/UnitCellVol;
  Vector b3 = (a1.Cross(a2))/UnitCellVol;

  cout << "UnitCellVol:\t" << UnitCellVol << endl;

  LatticePlane Plane( b1, b2, b3, 0, 0, 1, &IronFormFactor, UnitCellVol);
  LatticePlane PlaneB( b1, b2, b3, 0, 0, 2, &IronFormFactor, UnitCellVol);
  LatticePlane PlaneC( b1, b2, b3, 0, 1, 1, &IronFormFactor, UnitCellVol);
  LatticePlane PlaneD( b1, b2, b3, 1, 1, 1, &IronFormFactor, UnitCellVol);

  float Wavelength0 = a0*sqrt(2.0f); //45 degrees

  cout << "|H|/2 " << Plane.h << Plane.k << Plane.l << ":\t" << Plane.HalfMagH << endl;
  cout << "|H|/2 " << PlaneB.h << PlaneB.k << PlaneB.l << ":\t" << PlaneB.HalfMagH << endl;
  cout << "|H|/2 " << PlaneC.h << PlaneC.k << PlaneC.l << ":\t" << PlaneC.HalfMagH << endl;
  cout << "|H|/2 " << PlaneD.h << PlaneD.k << PlaneD.l << ":\t" << PlaneD.HalfMagH << endl;

  float ProgressCounter = 0.0f;

  //8A->53deg
  //2A->11.5deg
  
  Wavelength0 = 5.0f;

  float WavelengthDeltaMin = -4.0f;
  float WavelengthDeltaMax = +5.0f;

  float WavelengthDeltaStep = 0.05f;

  float fNumWavelengthPoints = ( (WavelengthDeltaMax - WavelengthDeltaMin)/WavelengthDeltaStep ) + 1.0f;

  float IntNumPoints = float(int(fNumWavelengthPoints)); //need to write out a float for gnuplot

  cout << "fNumWavelengthPoints:\t" << fNumWavelengthPoints << endl;
  cout << "IntNumPoints:\t" << IntNumPoints << endl;

  //ofstream ResultsFile;
  //ResultsFile.open( "Results.txt" );

  //ofstream ResultsBinary;
  //ResultsBinary.open ("ResultsBinary", ios::out | ios::binary);

  //ResultsBinary.write( (char*)&(IntNumPoints), sizeof(IntNumPoints) );

  //ResultsFile << IntNumPoints << "\t";

  /* for( int i = 0; i < IntNumPoints; i++)
  {
    float Wavelength = (Wavelength0 + WavelengthDeltaMin) + WavelengthDeltaStep*float(i);
    ResultsBinary.write( (char*)&(Wavelength), sizeof(Wavelength) );
    //ResultsFile << IntNumPoints << "\t";
    }*/
  
  


  float WidthLimit = Deg2Rad(3.0f);

  for( float WavelengthDelta = WavelengthDeltaMin; 
	 WavelengthDelta <= WavelengthDeltaMax;
	 WavelengthDelta += WavelengthDeltaStep)
  {

    float Wavelength = Wavelength0 + WavelengthDelta;

    float BraggA = Plane.FindBraggReflectionAngle( Wavelength );
    float BraggB = PlaneB.FindBraggReflectionAngle( Wavelength );
    float BraggC = PlaneC.FindBraggReflectionAngle( Wavelength );
    float BraggD = PlaneD.FindBraggReflectionAngle( Wavelength );

    //TODO: Check for bad bragg angles here?

    float WidthA = Plane.ScherrerWidth( Wavelength, 100);  WidthA = UpperLimit( WidthA, WidthLimit);    
    float WidthB = PlaneB.ScherrerWidth( Wavelength, 100); WidthB = UpperLimit( WidthB, WidthLimit);
    float WidthC = PlaneC.ScherrerWidth( Wavelength, 100); WidthC = UpperLimit( WidthC, WidthLimit);
    float WidthD = PlaneD.ScherrerWidth( Wavelength, 100); WidthD = UpperLimit( WidthD, WidthLimit);

    float CumulativeI = 0.0f;

    for( float Theta = 0.0f; Theta <= 90.0f; Theta += 0.001)    //dtheta 0.01
    {

    //ResultsBinary.write( (char*)&(Theta), sizeof(Theta) );

      //Do this in radians
   
      float ThetaRad = Deg2Rad(Theta);
      float I = Plane.CalculateIntensitySimple( Wavelength, ThetaRad, WidthA, BraggA );
      I += PlaneB.CalculateIntensitySimple( Wavelength, ThetaRad, WidthB, BraggB );
      I += PlaneC.CalculateIntensitySimple( Wavelength, ThetaRad, WidthC, BraggC );
      I += PlaneD.CalculateIntensitySimple( Wavelength, ThetaRad, WidthD, BraggD );

      CumulativeI += I;
      
      //if( I > 0.5f )
      //{
      //ResultsFile << Theta << '\t' << Wavelength << '\t' <<  CumulativeI << endl;	
	//}
      
      //ResultsBinary.write( (char*)&(I), sizeof(I) );
      

      //Integral += double(I) * double(0.0001);

      //ResultsFile << Theta << '\t' << Wavelength << '\t' <<  I << endl;
      //ResultsFile << Theta << '\t' << I << endl;  
    }

    //ResultsFile << endl;// For pm3d


    float Progress = (WavelengthDelta + 4.0f )/(9.0f);

    if(Progress >= ProgressCounter)
    {
      cout << "Progress:\t" << Progress * 100.0f << endl;
      ProgressCounter += 0.1f;
    }


  }

  //ResultsFile.close();

  //ResultsBinary.close();

  /*ifstream InputBinary;
  InputBinary.open ("ResultsBinary", ios::in | ios::binary);

  float NumPoints = 0.0;

  InputBinary.read( (char*)&NumPoints, sizeof(float) );

  cout << NumPoints << endl;

  InputBinary.close();*/

  return 0;
}
