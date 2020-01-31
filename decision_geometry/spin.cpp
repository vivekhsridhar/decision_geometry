//
//  spin.cpp
//  decision_geometry
//
//  Created by Vivek Hari Sridhar on 04/12/19.
//  Copyright © 2019 Vivek Hari Sridhar. All rights reserved.
//

#include "spin.h"

spin::spin(void)
{
}

spin::~spin(void)
{
}

void spin::Setup(const CVec2D& set_position, double& set_temperature, int& set_informed, bool& set_state, bool& set_picked)
{
    position = set_position;
    informed = set_informed;
    state = set_state;
    picked = set_picked;
    
    temperature = set_temperature;
}

void spin::AddPreference(CVec2D& cue_centre)
{
    preference = (cue_centre - position).normalise();
}

int spin::GetInformed()
{
    return informed;
}

void spin::SetInformed(int& info)
{
    informed = info;
}

void spin::Copy(spin& source)
{
    preference = source.preference;
    informed = source.informed;
    state = source.state;
    picked = source.picked;
    temperature = source.temperature;
}
