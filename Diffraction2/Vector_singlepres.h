#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <cmath>
#include <iostream>
#include <stdio.h>
using namespace std;


#define PI (3.14159265)

float Deg2Rad( float Degrees )
{
  return Degrees * (PI/180.0f);
}

float Rad2Deg( float Radians )
{
  return Radians * (180.0f/PI);
}


float Rad2Deg2( float Radians )
{
  if( Radians <= 0.0f)
  {
    return Radians;
  }
  else
  {
    return Radians * (180.0f/PI);
  }
}

float Deg2Rad2( float Degrees )
{
  if( Degrees <= 0.0f)
  {
    return Degrees;
  }
  else
  {
    return Degrees * (PI/180.0f);
  }
}



class Vector
{
public:
  float x,y,z;

  Vector ();
  Vector ( float X, float Y, float Z);
  Vector ( float r, float theta, float phi, bool bSpherical)
  {
    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);
  }

  // slightly more efficient way to create a spherical unit vector.
  //todo: get around compiler error
  //Need the extra input so that there's no compiler error.
  //(float,float,bool) is apparently equivalent to (float,float,float)
  //might be a GNU thing?
  Vector ( float theta, float phi, bool bSpherical, bool bUnit)
  {
    x = sin(theta)*cos(phi);
    y = sin(theta)*sin(phi);
    z = cos(theta);
  }

  float GetTheta()
  {
    return acos(z/Magnitude());
  }

  void SetTheta( float theta )
  {
    float r = Magnitude();
    float phi = GetPhi();

    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);
  }

  float GetPhi()
  {
    return atan2(y,x);
  }

  void SetPhi( float phi)
  {
    float r = Magnitude();
    float theta = GetTheta();

    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);
  }

  float Dot( const Vector V )
  {
    return x*V.x + y*V.y + z*V.z;
  }

  Vector Cross ( const Vector V)
  {
    return Vector( y*V.z - z*V.y , z*V.x - x*V.z , x*V.y - y*V.x);
  }

  float Magnitude()
  {
    return sqrt(x*x + y*y + z*z);
  }

  Vector Normalized()
  {
    float Mag = this->Magnitude();
    if(Mag == 0.0)
    {
      return Vector(0,0,0); //avoid the divide by 0 error
    }

    return this->Scale(float(1.0)/Mag);
  }

  Vector Scale( const float Scale )
  {
    Vector Out;
    Out.x = x*Scale;
    Out.y = y*Scale;
    Out.z = z*Scale;
    return Out;
  }

  Vector operator+ (Vector V)
  {
    return this->Add(V);
  }

  Vector operator- (Vector V)
  {
    return this->Subtract(V);
  }

  Vector operator* (float Scale)
  {
    return this->Scale( Scale );
  }


  Vector Divide( float Denominator )
  {
    Vector Out;
    Out.x = x/Denominator;
    Out.y = y/Denominator;
    Out.z = z/Denominator;
    return Out;
  }

  Vector operator/ (float Denominator)
  {
    return this->Divide( Denominator );
  }


  Vector Add( Vector V)
  {
    Vector Out;
    Out.x = x + V.x;
    Out.y = y + V.y;
    Out.z = z + V.z;
    return Out;
  }

  Vector Subtract( Vector V)
  {
    Vector Out;
    Out.x = x - V.x;
    Out.y = y - V.y;
    Out.z = z - V.z;
    return Out;
  }

  void Print(bool newline = true)
  {
    if(newline)
    {
      printf("%f\t%f\t%f\n", x,y,z);
    }
    else
    {
      printf("%f\t%f\t%f", x,y,z);
    }
  }

};

Vector operator*(float Scale, Vector V)
{
  return V*Scale;
}

Vector::Vector(float x, float y, float z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

Vector::Vector(void)
{
  x = 0;
  y = 0;
  z = 0;
}


//MUST be a unit vector. Uses euler matrix or from monte carlo paper (same but phi=0 in euler)
Vector ScatterUnitVector( Vector In, float Theta, float Phi)
{
  /*
    cos(theta)cos(phi)     -sin(phi)     sin(theta)cos(phi)
    cos(theta)sin(phi)     cos(phi)      sin(theta)sin(phi)
    -sin(theta)               0               cos(theta)
   */

  Vector Out;

  float CT = cos(Theta);
  float ST = sin(Theta);
  float CP = cos(Phi);
  float SP = sin(Phi);


  /*
  Out.x = CT*CP*In.x - SP*In.y + ST*CP*In.z;
  Out.y = CT*SP*In.x + CP*In.y + ST*SP*In.z;
  Out.z = CT*In.z - ST*In.x;*/


  //http://en.wikipedia.org/wiki/Monte_Carlo_method_for_photon_transport#Step_3:_Absorption_and_scattering
  //And: NUCLEAR SCIENCE AND ENGINEERING: 131, 132–136 ~1999 Direction Cosines and Polarization Vectors
  //for Monte Carlo Photon Scattering

  float zz = sqrt(1.0-In.z*In.z);

  if( fabs(zz) <= 0.001 )
  {
    Out.x = ST*CP;
    Out.y = ST*SP;
    Out.z = CT*In.z/fabs(In.z);
    return Out;
  }

  Out.x = In.x*CT + (1.0/zz)*ST*(In.x*In.z*CP-In.y*SP);
  Out.y = In.y*CT + (1.0/zz)*ST*(In.y*In.z*CP+In.x*SP);
  Out.z = In.z*CT - CP*ST*zz;

  return Out;



}


#endif
