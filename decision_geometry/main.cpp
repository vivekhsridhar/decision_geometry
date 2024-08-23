//
//  main.cpp
//  decision_geometry
//
//  Created by Vivek Sridhar on 04/12/19.
//  Copyright © 2019 Vivek Sridhar. All rights reserved.
//

#include "parameteres.h"

using namespace rnd;

std::ofstream outputFile1;
std::ofstream outputFile2;
//std::ofstream outputFile3;

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
    symmetric = false;
    lek = true;
    
    // time parameters
    num_replicates = 100;
    num_timesteps = 100000;
    timestep_number = 0;
    trial_time = 0;
    
    // space parameters
    arena_size = 600;
    if (number_of_cues == 2) max_angle = PI/36;
    else max_angle = 4*PI/9;
    
    start_dist = 500.0;
    dist_thresh = 2.0;
    
    if (!lek) arena_centre = CVec2D((double)arena_size / 2, (double)arena_size / 2);
    else arena_centre = CVec2D(500.0, 500.0);
    
    // system parameters
    total_agents = 240;
    n_flips = 0;
    nu = 0.5;
    dev = 0.02;
    system_energy = 0.0;
    system_magnetisation = CVec2D(0.0, 0.0);
    
    // run parameters
    reset_no = 0;
    std::fill_n(n_inds_preference, number_of_cues, 0);
    
    // output parameter
    energy = 0.0;
    magnetisation = CVec2D(0.0, 0.0);
    preference = new CVec2D[total_agents+1];
    state = new bool[total_agents+1];
    
    // output variables
    cue_reached = -1;
    output_frequency = 10;
    path_length = 0.0;
    centroid = CVec2D(0.0,0.0);
    filename = "/Users/vivekhsridhar/Library/Mobile Documents/com~apple~CloudDocs/Documents/Code/decision_geometry/decision_geometry/output/territories.csv";
    
    // class vectors
    agent = new spin[total_agents];
    CS = new cue[number_of_cues];
    
    // open output files
    outputFile1 = std::ofstream("geometry.csv");
    outputFile2 = std::ofstream("targets.csv");
    //outputFile3 = std::ofstream("target_reached.csv");
    
    // output file headers
    outputFile1 << "time" << ", " << "x" << ", " << "y" << ", " << "angular_difference" << ", " << "direction_chosen" << "\n";
    outputFile2 << "target_id" << ", " << "target_x" << ", " << "target_y" << "\n";
    //outputFile3 << "target_reached" << ", " << "trial_duration" << ", " << "speed" << "\n";
    
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
    double temp = 0.0;

    SetupSimulation(temp);
    for (int rep = 0; rep != num_replicates; ++rep)
    {
        ResetSetup(0, arena_centre.y);
        
        while (trial_time != num_timesteps)
        {
            FlipSpins();
            MoveAgents(temp);
            if (trial_time % output_frequency == 0 && trial_time != 0)
            {
                CalculateMagnetisation();
                GenerationOutput(rep);
            }
            
            ++trial_time;
            ++timestep_number;
            
            // reset agents if target is reached
            if (rep_done) break;
        }
        
        if (rep % 20 == 0) std::cout << rep << " ";
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
        agent[i].position += system_magnetisation;
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
            
            //outputFile3 << i << ", " << trial_time << ", " << path_length / trial_time << "\n";
        }
    }
    
    path_length += system_magnetisation.length();
}

void CalculateEnergy(int spin_id)
{
    energy = 0.0;
    for (int i = 0; i != total_agents+1; ++i)
    {
        double ang = preference[spin_id].smallestAngleTo(preference[i]) * PiOver180;
        
        ang = PI * pow(ang / PI, nu);
        double J = cos(ang);
        
        if (i != spin_id) energy -=  J * state[spin_id] * state[i];
    }
    energy /= (total_agents+1);
}

void CalculateMagnetisation()
{
    magnetisation = CVec2D(0.0, 0.0);
    for (int i = 0; i != total_agents+1; ++i) magnetisation += preference[i] * state[i];
    magnetisation /= (total_agents+1);
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
//    if (lek) SetupEnvironmentRandom();
    if (lek) SetupEnvironmentFromCSV(filename);
    else if (symmetric) SetupEnvironmentSymmetric();
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

void SetupEnvironmentFromCSV(const std::string& file_name)
{
    std::vector<CVec2D> positions;

    // Open the CSV file
    std::ifstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_name << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); // Read and ignore the header line

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string idx_str, x_str, y_str;

        // Read the idx, pos_x, and pos_y from the line
        std::getline(ss, idx_str, ','); // Read idx
        std::getline(ss, x_str, ',');   // Read pos_x
        std::getline(ss, y_str, ',');   // Read pos_y

        // Convert pos_x and pos_y to double
        double x = std::stod(x_str);
        double y = std::stod(y_str);

        // Create a CVec2D object and add it to the vector
        positions.emplace_back(x, y);
    }

    file.close();

    // Verify if we have the correct number of positions
    if (positions.size() != number_of_cues) {
        std::cerr << "Error: Number of positions read (" << positions.size() << ") does not match the number of cues (" << number_of_cues << ")." << std::endl;
        return;
    }

    // Set up the environment with the read positions
    for (std::size_t i = 0; i < positions.size(); ++i) {
        centres[i] = positions[i];
        CS[i].Setup(centres[i]);
        
        outputFile2 << i << ", " << centres[i].x << ", " << centres[i].y << "\n";
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
        ++n_inds_preference[set_informed];
        
        set_temperature = temp;
        set_deviation = rnd::normal(0.0, dev);
            
        agent[i].Setup(set_position, set_temperature, set_informed, set_state, set_deviation, set_picked);
        agent[i].AddPreference(CS[agent[i].GetInformed()].centre);
        
        preference[i] = agent[i].preference;
        state[i] = agent[i].state;
    }
    preference[total_agents] = (CS[0].centre - centroid).normalise();
    state[total_agents] = true;
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
        agent[i].prime_deviation = rnd::normal(0.0, dev);
        
        preference[i] = agent[i].preference;
        state[i] = true;
        
        preference[total_agents] = (CS[rnd::integer(number_of_cues)].centre - centroid).normalise();
        state[total_agents] = true;
    }
    
    trial_time = 0;
    path_length = 0.0;
    n_flips = 0;
    
    ++reset_no;
    
    rep_done = false;
}

void ResetStates()
{
    for (int i = 0; i != total_agents; ++i)
    {
        preference[i] = agent[i].preference;
        state[i] = agent[i].state;
    }
    preference[total_agents] = (CS[0].centre - centroid).normalise();
    state[total_agents] = true;
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
    CVec2D v0;
    CVec2D v1;
    CVec2D v2;
    v0 = system_magnetisation.normalise();
    v1 = (CS[0].centre - centroid).normalise();
    v2 = (CS[number_of_cues - 1].centre - centroid).normalise();
    
    
    outputFile1 << trial_time << ", " << centroid.x << ", " << centroid.y << ", " << v1.smallestAngleTo(v2) << ", " << v0.smallestAngleTo(v1) << "\n";
}
