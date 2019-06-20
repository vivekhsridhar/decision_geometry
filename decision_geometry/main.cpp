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
std::ofstream outputFile3;
std::ofstream outputFile4;

//**************************************************************************************************
//**	MAIN	************************************************************************************
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
    arena_size = 1000;
    arena_centre = CVec2D((double)arena_size / 2, (double)arena_size / 2);
    total_agents = 6;
    system_energy = 0.0;
    system_magnetisation = CVec2D(0.0, 0.0);
    
    start_dist = 500.0;
    overall_angle = 360.0;  // used for the symmetric case ('overall_angle' is split into 'number_of_cues' equal angles)
    dist_thresh = 10.0;
    max_angle = Pi/3;       // used for the asymmetric case ('max_angle' is split into 'number_of_cues' - 1 equal angles)
    
    agent = new spin[total_agents];
    CS = new cue[number_of_cues];
    
    std::fill_n(n_inds_preference, number_of_cues, 0);
    symmetric = false;
    
    // Open output files
    outputFile1 = std::ofstream("geometry.csv");
    outputFile2 = std::ofstream("success_rate.csv");
    outputFile3 = std::ofstream("dtime.csv");
    outputFile4 = std::ofstream("speed_profile.csv");
    
    // Output file headers
    outputFile1 << "temperature" << ", " << "replicate" << ", " << "x" << ", " << "y" << ", " << "theta" << ", " << "dir_x" << ", " << "dir_y" << "\n";
    outputFile2 << "temperature" << ", " << "replicate" << ", " << "ncues" << ", " << "n_success" << "\n";
    outputFile3 << "group_size" << ", " << "ncues" << ", " << "path_length" << ", " << "decision_time" << "\n";
    outputFile4 << "time" << ", " << "speed" << "\n";
    
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
    int num_timesteps = 100000;
    int num_replicates = 1;
    timestep_number = 0;
    
    for (double temp = 0.0; temp < 0.1;)
    {
        SetupSimulation(temp);
        for (int rep = 0; rep != num_replicates; ++rep)
        {
            ResetSetup();
            while (timestep_number != num_timesteps)
            {
                FlipSpins();
                MoveAgents();
                if (timestep_number % 10 == 0)
                {
                    Graphics();
                    GenerationalOutput(temp, rep);
                }
                
                ++trial_time;
                ++timestep_number;
                
                // reset agents if target is reached
                if (rep_done)
                {
                    ++n_successful;
                    break;
                }
            }
            
            outputFile2 << temp << ", " << rep << ", " << number_of_cues << ", " << n_successful << "\n";
        }
        
        temp += 0.02;
        std::cout << temp << " ";
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
        double J = agent[spin_id].preference.dot(agent[i].preference);
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

void MoveAgents()
{
    for (int i = 0; i != total_agents; ++i)
    {
        agent[i].position += system_magnetisation;
        agent[i].AddPreference(CS[agent[i].GetInformed()].centre, arena_centre);
    }
    
    for (int i = 0; i != number_of_cues; ++i)
    {
        if (centroid.distanceTo(CS[i].centre) < dist_thresh * dist_thresh)
        {
            rep_done = true;
            cue_reached = i;
            for (int j = 0; j != total_agents; ++j)
            {
                if (agent[j].GetInformed() == i) agent[j].fitness += 1.0 / n_inds_preference[i];
            }
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
    n_successful = 0;
    trial_time = 0;         // time elapsed within trial
    reset_no = 0;           // trial number
    cue_reached = -1;
    
    avg_speed = 0.0;        // average speed of the group through the trial
    path_length = 0.0;      // length of the path the group takes (a proxy for how convoluted the path is)
    
    centroid = arena_centre;
    if (symmetric) SetupEnvironmentSymmetric();
    else SetupEnvironmentAsymmetric();
    
    SetupSpins(temp);
}

void SetupEnvironmentSymmetric()
{
    CVec2D start;
    start = arena_centre;
    
    double theta = overall_angle * PiOver180 / number_of_cues;
    
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

void SetupSpins(double temp)
{
    CVec2D set_position;
    CVec2D set_preference;
    int set_informed;
    bool set_state;
        
    double set_temperature;
    double set_fitness;
    
    set_preference = CVec2D(0.0, 0.0);
    set_fitness = 0.0;
    
    for(int i = 0; i != total_agents; ++i)
    {
        set_position = RandomBoundedPoint();
            
        if (rnd::uniform() < 0.0) set_state = false;
        else set_state = true;
            
        set_informed = i % number_of_cues;
        ++n_inds_preference[set_informed];
            
        set_temperature = temp;
            
        agent[i].Setup(set_position, set_temperature, set_informed, set_state, set_fitness);
        agent[i].AddPreference(CS[agent[i].GetInformed()].centre, arena_centre);
    }
}

void ResetSetup()
{
    outputFile3 << total_agents << ", " << number_of_cues << ", " << path_length << ", " << trial_time << "\n";
    
    for(int i = 0; i != total_agents; ++i)
    {
        agent[i].position = RandomBoundedPoint();
        
        if (rnd::uniform() < 0.5) agent[i].state = false;
        else agent[i].state = true;
        
        int info = i % number_of_cues;
        agent[i].SetInformed(info);
        agent[i].preference = CVec2D(0.0, 0.0);
    }
    
    if (symmetric) SetupEnvironmentSymmetric();
    else SetupEnvironmentAsymmetric();
    
    trial_time = 0;
    path_length = 0.0;
    avg_speed = 0.0;
    
    ++reset_no;
    
    rep_done = false;
}

CVec2D RandomBoundedPoint()
{
    // Create randomly distributed co-ordinates in the simulated world space
    double range_x = (double) arena_size;
    double range_y = (double) arena_size;
    
    double random_x = uniform();
    double random_y = uniform();
    
    // Individuals start in the centre-left 100th of their world
    random_x *= (range_x / 100.0);
    random_y *= (range_y / 100.0);
    if (symmetric) random_x += (double) (range_x / 2.0 - range_x / 200.0);
    random_y += (double) (range_y / 2.0 - range_y / 200.0);
    CVec2D random_point(random_x, random_y);
    return random_point;
}

//**************************************************************************************************
//**    OUTPUT  ************************************************************************************
//**************************************************************************************************

void GenerationalOutput(double temp, int rep)
{
    CVec2D v1;
    CVec2D v2;
    v1 = (CS[0].centre - centroid).normalise();
    v2 = (CS[number_of_cues-1].centre - centroid).normalise();
    
    outputFile1 << temp << ", " << rep << ", " << centroid.x << ", " << centroid.y << ", " << v1.smallestAngleTo(v2) << ", " << system_magnetisation.x << ", " << system_magnetisation.y << "\n";
    outputFile4 << trial_time << ", " << system_magnetisation.length() << "\n";
}

//**************************************************************************************************
//**	GRAPHICS   *********************************************************************************
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

void GraphicsWriter(VideoWriter& video_writer, int& timestep_number, const int& num_timesteps)
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
        int colour = agent[i].state;
        circle(visualisation, Point(agent[i].position.x, agent[i].position.y), 2, colours[colour], -1);
    }
    
    // Display timestep number & cue counter on screen
    putText(visualisation, std::to_string(timestep_number), Point(10, 30), FONT_HERSHEY_COMPLEX_SMALL, 0.8, Scalar(200, 200, 250), 1, CV_AA);
    
    // Write video
    video_writer.write(visualisation); //writer the frame into the file
    
    if (timestep_number == num_timesteps)
    {
        video_writer.release();
    }
    waitKey(1);
}
