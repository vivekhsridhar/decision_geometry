//
//  vector2D.h
//  decision_geometry
//
//  Created by Vivek Hari Sridhar on 04/12/19.
//  Copyright © 2019 Vivek Hari Sridhar. All rights reserved.
//

#ifndef vector2D_h
#define vector2D_h


#include <cmath>
#include <stdlib.h>
#include <sstream>
#include "random.h"

const double PiOver180 = 1.74532925199433E-002;
const double PiUnder180 = 5.72957795130823E+001;

class CVec2D
{
public:
    CVec2D(void);
    ~CVec2D(void);
    
    CVec2D(double x1, double x2);
    CVec2D(const CVec2D&);
    CVec2D(CVec2D&&) noexcept;
    double x, y;
    
    CVec2D operator+(const CVec2D&);
    CVec2D operator+=(const CVec2D&);
    CVec2D operator-(const CVec2D&);
    CVec2D operator-=(const CVec2D&);
    CVec2D operator-();
    CVec2D operator=(const CVec2D&);
    CVec2D operator*(const CVec2D&);
    CVec2D operator*=(const CVec2D&);
    CVec2D operator*(double);
    CVec2D operator*=(double);
    CVec2D operator/(const CVec2D&);
    CVec2D operator/(double);
    CVec2D operator/=(double);
    
    CVec2D normalise();
    double dot(const CVec2D&);
    double cross(const CVec2D&);
    void rotate(double degrees);
    double length();
    double distanceTo(const CVec2D&);
    inline void Clear() { x = 0.0; y = 0.0; }
    double polarAngle();
    double polarAngleZeroNorth();
    double smallestAngleTo(CVec2D);
};

#endif /* vector2D_h */
