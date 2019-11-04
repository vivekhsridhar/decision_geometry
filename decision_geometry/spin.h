//
//  spin.h
//  multi-choice_decision_geometry
//
//  Created by Vivek Hari Sridhar on 18/09/17.
//  Copyright © 2017 Vivek Hari Sridhar. All rights reserved.
//

#ifndef spin_h
#define spin_h
#include <iostream>
#include "cue.h"

class spin
{
public:
    spin(void);
    ~spin(void);
    
    void Setup(const CVec2D& set_position, double& set_temperature, int& set_informed, bool& set_state);
    void AddPreference(CVec2D& cue_centre);
    int GetInformed();
    void SetInformed(int& informed);
    void Copy(spin& source);
    
    CVec2D position;
    CVec2D preference;
    bool state;
    
    double temperature;
private:
    int informed;
};

#endif /* spin_h */
