//
//  main.cpp
//  multi-choice_decision_geometry
//
//  Created by Vivek Sridhar on 18/09/17.
//  Copyright © 2017 Vivek Sridhar. All rights reserved.
//

#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/videoio.hpp"
#include "parameteres.h"

using namespace rnd;
using namespace cv;

std::ofstream outputFile1;
std::ofstream outputFile2;

//**************************************************************************************************
//**    MAIN    ************************************************************************************
//**************************************************************************************************

int main()
{
    echo("simulation started");
    
    //===================================
    //==    model parameterisation  =====
    //===================================
    // Random generator engine from a time-based seed
    unsigned seed = static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    rnd::set_seed(seed);
    
    // Set model parameters
    distance = false;
    symmetric = false;
    assert(symmetric == false || distance == false);
    
    arena_size = 1000;
    arena_centre = CVec2D((double)arena_size / 2, (double)arena_size / 2);
    total_agents = 60;
    system_energy = 0.0;
    system_magnetisation = CVec2D(0.0, 0.0);
    
    left_right_dist = 500;
    if (distance) start_dist = 100.0;
    else start_dist = 500;
    
    dist_thresh = 10.0;
    if (number_of_cues == 2) max_angle = Pi/3;
    else max_angle = 4*Pi/9;    // used for the asymmetric case ('max_angle' is split into 'number_of_cues' - 1 equal angles)
    
    hat_width = 0.25;           // inversely proportional to the hat width in the mexican hat function
    path_length = 0.0;
    
    agent = new spin[total_agents];
    CS = new cue[number_of_cues];
    
    std::fill_n(n_inds_preference, number_of_cues, 0);
    
    // Open output files
    outputFile1 = std::ofstream("geometry.csv");
    outputFile2 = std::ofstream("cue_reached.csv");
    
    // Output file headers
    if (distance) outputFile1 << "x" << ", " << "y" << ", " << "left_right_distance" << ", " << "front_back_distance" << "\n";
    else outputFile1 << "replicate" << ", " << "x" << ", " << "y" << ", " << "angular_disagreement" << ", " << "dir_x" << ", " << "dir_y" << "\n";
    
    outputFile2 << "replicate" << ", " << "n1" << ", " << "n2" << ", ";
    if (number_of_cues == 3) outputFile2 << "n3" << ", ";
    outputFile2 << "target_reached" << ", " << "path_length" << "\n";
    
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
    int num_timesteps = 10000;
    int num_replicates = 50;
    timestep_number = 0;
    
    for (left_right_dist = 50; left_right_dist < 500; )
    {
        SetupSimulation(0.05);
        for (int rep = 0; rep != num_replicates; ++rep)
        {
            ResetSetup(0, arena_size/2);
            
            while (trial_time != num_timesteps)
            {
                FlipSpins();
                MoveAgents(rep);
                if (trial_time != 0 && trial_time % 10 == 0)
                {
                    //Graphics();
                    GenerationalOutput(rep);
                }
                
                ++trial_time;
                ++timestep_number;
                
                // reset agents if target is reached
                if (rep_done) break;
            }
        }
        
        if (distance) left_right_dist += 50;
        else left_right_dist += 500;
        std::cout << left_right_dist << " ";
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
        double ang = agent[spin_id].preference.smallestAngleTo(agent[i].preference) * PI / 180.0;
        double J = 1.8*(1 - hat_width * ang * ang) * exp(-hat_width * ang * ang) - 1.0;
        
        if (i != spin_id) system_energy -=  J * agent[spin_id].state * agent[i].state;
    }
    system_energy /= total_agents;
    
    // calculate magnetisation
    centroid = CVec2D(0.0, 0.0);
    system_magnetisation = CVec2D(0.0, 0.0);
    for (int i = 0; i != total_agents; ++i)
    {
        centroid += agent[i].position;
        system_magnetisation += agent[i].preference * agent[i].state;
    }
    centroid /= total_agents;
    system_magnetisation /= total_agents;
}

void MoveAgents(int rep)
{
    for (int i = 0; i != total_agents; ++i)
    {
        if (!distance) agent[i].position += system_magnetisation;
        else agent[i].position.y += system_magnetisation.y;
        agent[i].AddPreference(CS[agent[i].GetInformed()].centre);
    }
    
    for (int i = 0; i != number_of_cues; ++i)
    {
        if (centroid.distanceTo(CS[i].centre) < dist_thresh * dist_thresh)
        {
            rep_done = true;
            cue_reached = i;
            
            outputFile2 << rep << ", ";
            for (int j = 0; j != number_of_cues; ++j)
            {
                outputFile2 << n_inds_preference[j] << ", ";
            }
            outputFile2 << i << ", " << path_length << "\n";
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
    trial_time = 0;         // time elapsed within trial
    reset_no = 0;           // trial number
    cue_reached = -1;
    path_length = 0.0;      // length of the path the group takes (a proxy for how convoluted the path is)
    
    centroid = arena_centre;
    if (symmetric) SetupEnvironmentSymmetric();
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

void SetupSpins(double set_temperature)
{
    CVec2D set_position;
    CVec2D set_preference;
    int set_informed;
    bool set_state;
        
    double set_fitness;
    
    set_preference = CVec2D(0.0, 0.0);
    set_fitness = 0.0;
    
    for(int i = 0; i != total_agents; ++i)
    {
        set_position = RandomBoundedPoint(0, 0);
            
        if (rnd::uniform() < 0.5) set_state = false;
        else set_state = true;
            
        set_informed = i % number_of_cues;
        ++n_inds_preference[set_informed];
            
        agent[i].Setup(set_position, set_temperature, set_informed, set_state);
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
    
    if (symmetric) SetupEnvironmentSymmetric();
    else if (distance) SetupEnvironmentDistances();
    else SetupEnvironmentAsymmetric();
    
    trial_time = 0;
    path_length = 0.0;
    
    ++reset_no;
    
    rep_done = false;
}

CVec2D RandomBoundedPoint(double x, double y)
{
    double random_x = x + uniform() - 0.5;
    double random_y = y + uniform() - 0.5;
    
    CVec2D random_point(random_x, random_y);
    if (symmetric) random_point = arena_centre;
    
    return random_point;
}

//**************************************************************************************************
//**    OUTPUT  ************************************************************************************
//**************************************************************************************************

void GenerationalOutput(int rep)
{
    CVec2D v1;
    CVec2D v2;
    v1 = (CS[0].centre - centroid).normalise();
    v2 = (CS[number_of_cues-1].centre - centroid).normalise();
    
    if (distance)
    {
        outputFile1 << centroid.x << ", " << centroid.y << ", " << left_right_dist << ", " << start_dist << "\n";
    }
    else
    {
        outputFile1 << rep << ", " << centroid.x << ", " << centroid.y << ", " << v1.smallestAngleTo(v2) << ", " << system_magnetisation.x << ", " << system_magnetisation.y << "\n";
    }
}

//**************************************************************************************************
//**    GRAPHICS   *********************************************************************************
//**************************************************************************************************

void Graphics()
{
    // Colours vector
    Scalar colours[6] = {Scalar(234, 174, 84), Scalar(24, 202, 247), Scalar(60, 76, 231), Scalar(113, 204, 46), Scalar(34, 126, 230), Scalar(241, 240, 236)};
    
    // Draw area
    Mat visualisation = Mat::zeros(arena_size, arena_size, CV_8UC3);
    for (int i = 0; i != number_of_cues; ++i)
    {
        circle(visualisation, Point(CS[i].centre.x, CS[i].centre.y), 8, colours[i], -1, CV_AA);
    }
    
    // Draw spins
    for (int i = 0; i != total_agents; ++i)
    {
        int colour = agent[i].GetInformed();
        circle(visualisation, Point(agent[i].position.x, agent[i].position.y), 2, colours[colour], -1);
    }
    
    // Display timestep number & cue counter on screen
    putText(visualisation, std::to_string(timestep_number), cvPoint(10,30), FONT_HERSHEY_COMPLEX_SMALL, 0.8, cvScalar(200,200,250), 1, CV_AA);
    
    imshow("ising_model", visualisation);
    waitKey(1);
}
