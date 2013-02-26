#include <iostream>
#include "CCD.h"
#include "Vector.h"
#include <cstdlib>

using namespace std;

float UniformRand()
{
    return float(rand()) / (float(RAND_MAX)+1.0);
}

int main()
{
    Vector InputCCDOrigin(0,0,1); //origin of CCD
    Vector InputCCDNormal(0.0,0,1); //direction that CCD points in.
    double InputCCDAngle = 0;

    CCD CCDCamera(InputCCDOrigin, InputCCDNormal, InputCCDAngle,
                  1.0, 1.0);

    double RayLength;
    Vector IntersectPoint;
    double XIntersect, YIntersect;
    int XPixel,YPixel;

    Vector RayDirection(1,0,1);
    Vector RaySource(0,0,0);

    RayDirection = RayDirection.Normalized();

    if(CCDCamera.GetPixelRayCCDIntersect(RaySource,RayDirection,
                                         RayLength,IntersectPoint,
                                         XIntersect,YIntersect,
                                         XPixel,YPixel) == false)
    {
        cout << "No Intersection!" << endl;
        exit(0);
    }

    cout << "Real Intersect Point:\t";
    IntersectPoint.Print();
    cout << "Ray path length:\t" << RayLength << endl;
    cout << XIntersect << "\t" << YIntersect << endl;
    cout << XPixel << "\t" << YPixel << endl;

    return 0;
}
