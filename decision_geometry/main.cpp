//
//  main.cpp
//  decision_geometry
//
//  Created by Vivek Sridhar on 04/12/19.
//  Copyright © 2019 Vivek Sridhar. All rights reserved.
//

#include <cassert>
#include "parameteres.h"

using namespace rnd;

std::ofstream outputFile1;
std::ofstream outputFile2;
std::ofstream outputFile3;

//**************************************************************************************************
//**    MAIN    ************************************************************************************
//**************************************************************************************************

int main()
{
    echo("simulation started");
    
    // random generator engine from a time-based seed
    unsigned seed = static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    rnd::set_seed(seed);
    
    // boolean switches
    rep_done = false;
    symmetric = true;
    distance = false;
    lek = false;
    assert(distance == false || symmetric == false);
    assert(distance == false || lek == false);
    assert(symmetric == false || lek == false);
    
    // time parameters
    if (distance)
    {
        num_replicates = 50;
        num_timesteps = 10000;
    }
    else
    {
        num_replicates = 500;
        num_timesteps = 100000;
    }
    timestep_number = 0;
    trial_time = 0;
    equilibration_time = 1;
    
    // space parameters
    arena_size = 1000;
    if (number_of_cues == 2) max_angle = PI/3;
    else max_angle = 4*PI/9;
    
    if (distance) start_dist = 75.0;
    else start_dist = 500.0;
    dist_thresh = 10.0;
    left_right_dist = 500.0;
    centre_left_dist = 450;
    centre_right_dist = 150;
    arena_centre = CVec2D((double)arena_size / 2, (double)arena_size / 2);
    
    // system parameters
    total_agents = 60;
    n_flips = 0;
    nu = 0.5;
    A = 1.8;
    h = 0.25;
    c = 1.0;
    dev = 0.02;
    system_energy = 0.0;
    system_magnetisation = CVec2D(0.0, 0.0);
    assert(total_agents >= 5 * number_of_cues);
    
    // run parameters
    reset_no = 0;
    
    // output variables
    cue_reached = -1;
    output_frequency = 10;
    path_length = 0.0;
    centroid = CVec2D(0.0,0.0);
    
    // class vectors
    agent = new spin[total_agents];
    CS = new cue[number_of_cues];
    
    // open output files
    outputFile1 = std::ofstream("geometry.csv");
    if (!distance) 
    {
        outputFile2 = std::ofstream("targets.csv");
        outputFile3 = std::ofstream("targets_reached.csv");
    }
    
    // output file headers
    if (distance) outputFile1 << "time" << ", " << "x" << ", " << "y" << ", " << "left_right_distance" << ", " << "front_back_distance" << "\n";
    else
    {
        outputFile1 << "time" << ", " << "x" << ", " << "y" << "\n";
        outputFile2 << "target_id" << ", " << "target_x" << ", " << "target_y" << "\n";
        outputFile3 << "replicate" << ", " << "angle" << ", " << "target_reached" << ", " << "n_flips" << ", " << "trial_duration" << ", " << "speed" << "\n";
    }
    
    //===================================
    //==    functions in the main   =====
    //===================================
    RunGeneration();
    
    echo("simulation ended");
    return 0;
}

//**************************************************************************************************
//**    WITHIN GENERATIONAL FUNCTIONS   ************************************************************
//**************************************************************************************************

void RunGeneration()
{
    double temp = 0.1;
    for (left_right_dist = 0; left_right_dist <= 500; )
    {
        SetupSimulation(temp);
           
        for (int rep = 0; rep != num_replicates; ++rep)
        {
            ResetSetup(0, arena_size/2);
            
            while (trial_time != num_timesteps)
            {
                FlipSpins();
                MoveAgents(temp);
                if (trial_time % output_frequency == 0 && trial_time != 0) GenerationOutput(rep);
                
                ++trial_time;
                ++timestep_number;
                
                // reset agents if target is reached
                if (rep_done) break;
            }
            
            if (!distance && rep % 50 == 0) std::cout << rep << " ";
        }
        
        if (distance)
        {
            std::cout << left_right_dist << " " << start_dist << "\n";
            
            left_right_dist += 10;
            start_dist += 1.5;
        }
        else left_right_dist += 1000;
    }
}

void FlipSpins()
{
    int id = rnd::integer(total_agents);
    CalculateSystemProperties(id);
    double before = system_energy;
    agent[id].state = !agent[id].state;
    CalculateSystemProperties(id);
    double after = system_energy;
    
    double p_accept = 0.0;
    if (before < after) p_accept = exp(-(after - before) / agent[id].temperature);
    else p_accept = 1.0;
    
    if (rnd::uniform() >= p_accept) agent[id].state = !agent[id].state;
}

void CalculateSystemProperties(int spin_id)
{
    // calculate energy
    system_energy = 0.0;
    for (int i = 0; i != total_agents; ++i)
    {
        double ang = agent[spin_id].preference.smallestAngleTo(agent[i].preference) * PiOver180;
        
        ang = PI * pow(ang / PI, nu);
        double J = cos(ang);
//        double J = A * (1 - h * ang * ang) * exp(-h * ang * ang) - c;
        
        if (i != spin_id) system_energy -=  J * agent[spin_id].state * agent[i].state * agent[spin_id].picked * agent[i].picked;
    }
    system_energy *= number_of_cues;
    system_energy /= total_agents;
    
    // calculate magnetisation
    int total_picked = 0;
    centroid = CVec2D(0.0, 0.0);
    system_magnetisation = CVec2D(0.0, 0.0);
    for (int i = 0; i != total_agents; ++i)
    {
        centroid += agent[i].position;
        system_magnetisation += agent[i].preference * agent[i].state;
        total_picked += agent[i].picked;
    }
    centroid /= total_agents;
    system_magnetisation /= total_picked;
}

void MoveAgents(int rep)
{
    for (int i = 0; i != total_agents; ++i)
    {
        if (distance) agent[i].position.y += system_magnetisation.y;
        else agent[i].position += system_magnetisation;
        agent[i].AddPreference(CS[agent[i].GetInformed()].centre);
        
        double summation = 0.0;
        for (int j = 0; j != number_of_cues; ++j)
        {
            agent[i].GetDeviation(CS[j].centre, j);
            agent[i].probabilities[j] = GetProbability(agent[i].deviations[j], 0.0, dev);
            summation += agent[i].probabilities[j];
        }
        
        for (int j = 0; j != number_of_cues; ++j) agent[i].probabilities[j] /= summation;
        
        if (rnd::uniform() < agent[i].probabilities[agent[i].GetInformed()]) agent[i].picked = true;
        else agent[i].picked = false;
    }
    
    for (int i = 0; i != number_of_cues; ++i)
    {
        if (centroid.distanceTo(CS[i].centre) < dist_thresh * dist_thresh)
        {
            rep_done = true;
            cue_reached = i;
            
            outputFile3 << rep << ", " << max_angle << ", " << i << ", " << n_flips << ", " << trial_time << ", " << path_length / trial_time << "\n";
        }
    }
    
    path_length += system_magnetisation.length();
}

//**************************************************************************************************
//**    SETUP FUNCTIONS ****************************************************************************
//**************************************************************************************************

void SetupSimulation(double temp)
{
    timestep_number = 0;
    trial_time = 0;
    reset_no = 0;
    cue_reached = -1;
    
    path_length = 0.0;
    
    centroid = arena_centre;
    if (lek) SetupEnvironmentRandom();
    else if (symmetric) SetupEnvironmentSymmetric();
    else if (distance) SetupEnvironmentDistances();
    else SetupEnvironmentAsymmetric();
    
    SetupSpins(temp);
}

void SetupEnvironmentSymmetric()
{
    CVec2D start;
    start = arena_centre;
    
    double theta = 360.0 * PiOver180 / number_of_cues;
    
    for (int i = 0; i != number_of_cues; ++i)
    {
        centres[i] = start + CVec2D(start_dist * cos((i-1) * theta), start_dist * sin((i-1) * theta));
        CS[i].Setup(centres[i]);
        
        outputFile2 << i << ", " << centres[i].x << ", " << centres[i].y << "\n";
    }
}

void SetupEnvironmentAsymmetric()
{
    CVec2D start;
    start = CVec2D(0.0, arena_size / 2);
    
    double theta = 0.0;
    if (number_of_cues != 1) theta = max_angle / (number_of_cues - 1);
    
    for (int i = 0; i != number_of_cues; ++i)
    {
        centres[i] = start + CVec2D(start_dist * cos(i * theta - max_angle/2), start_dist * sin(i * theta - max_angle/2));
        CS[i].Setup(centres[i]);
        
        outputFile2 << i << ", " << centres[i].x << ", " << centres[i].y << "\n";
    }
}

void SetupEnvironmentDistances()
{
    for (int i = 0; i != number_of_cues; ++i)
    {
        if (number_of_cues == 2)
        {
            centres[i] = CVec2D(start_dist, arena_size / 2 - left_right_dist / 2 + i * left_right_dist);
        }
        else
        {
            centres[i] = CVec2D(start_dist, arena_size / 2 - left_right_dist + i * left_right_dist);
        }
        
        CS[i].Setup(centres[i]);
    }
}

void SetupEnvironmentRandom()
{
    for (int i = 0; i != number_of_cues; ++i)
    {
        centres[i] = arena_centre + RandomPolarPoint(0, arena_size/2 * arena_size/2);
        CS[i].Setup(centres[i]);
        
        outputFile2 << i << ", " << centres[i].x << ", " << centres[i].y << "\n";
    }
}

void SetupSpins(double temp)
{
    CVec2D set_position;
    CVec2D set_preference;
    int set_informed;
    bool set_state;
    bool set_picked;
    
    double set_temperature;
    double set_deviation;
    
    n_flips = 0;
    set_preference = CVec2D(0.0, 0.0);
    
    for(int i = 0; i != total_agents; ++i)
    {
        set_position = RandomBoundedPoint(0, 0);
            
        if (rnd::uniform() < 0.5) set_state = false;
        else set_state = true;
        set_picked = true;
        
        set_informed = i % number_of_cues;
        
        set_temperature = temp;
        set_deviation = rnd::normal(0.0, dev);
            
        agent[i].Setup(set_position, set_temperature, set_informed, set_state, set_deviation, set_picked);
        agent[i].AddPreference(CS[agent[i].GetInformed()].centre);
    }
}

void ResetSetup(double x, double y)
{
    for(int i = 0; i != total_agents; ++i)
    {
        agent[i].position = RandomBoundedPoint(x, y);
        
        if (rnd::uniform() < 0.5) agent[i].state = false;
        else agent[i].state = true;
        
        int info = i % number_of_cues;
        agent[i].SetInformed(info);
        agent[i].preference = CVec2D(0.0, 0.0);
    }
    
    trial_time = 0;
    path_length = 0.0;
    n_flips = 0;
    
    ++reset_no;
    
    rep_done = false;
}

CVec2D RandomBoundedPoint(double x, double y)
{
    double random_x;
    double random_y;
    CVec2D random_point;
    
    random_x = x + uniform() - 0.5;
    random_y = y + uniform() - 0.5;
    
    random_point = CVec2D(random_x, random_y);
    if (symmetric) random_point = arena_centre;
    
    return random_point;
}

CVec2D RandomPolarPoint(double rmin, double rmax)
{
    double random_r = (rmax - rmin) * uniform() + rmin;
    double random_theta = 2 * PI * uniform() - PI;
    double sqrt_r = sqrt(random_r);  // Avoid recalculating sqrt

    return CVec2D(sqrt_r * cos(random_theta), sqrt_r * sin(random_theta));
}

double GetProbability(double x, double mu, double sigma)
{
    double coeff = 1 / (sigma * sqrt(2*PI));
    double exp_function = exp(-pow(x - mu, 2) / (2 * pow(sigma, 2)));
    
    double prob = coeff * exp_function;
    return prob;
}

//**************************************************************************************************
//**    OUTPUT  ************************************************************************************
//**************************************************************************************************

void GenerationOutput(int rep)
{
    if (distance) outputFile1 << trial_time << ", " << centroid.x << ", " << centroid.y << ", " << left_right_dist << ", " << start_dist << "\n";
    else outputFile1 << trial_time << ", " << centroid.x << ", " << centroid.y << "\n";
}
