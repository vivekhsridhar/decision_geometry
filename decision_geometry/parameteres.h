//
//  parameteres.h
//  decision_geometry
//
//  Created by Vivek Hari Sridhar on 04/12/19.
//  Copyright © 2019 Vivek Hari Sridhar. All rights reserved.
//

#ifndef parameteres_h
#define parameteres_h

#include "spin.h"
#include <fstream>

// boolean switches
bool    rep_done;
bool    symmetric;
bool    lek;

// time parameters
int     num_replicates;
int     num_timesteps;
int     timestep_number;
int     trial_time;

// space parameters
int     arena_size;
double  max_angle;
double  start_dist;
double  dist_thresh;
double  overlap_sd;
CVec2D  arena_centre;

// system parameters
int     total_agents;
int     n_flips;
double  nu;
double  dev;
double  system_energy;
CVec2D  system_magnetisation;

// run parameters
int     reset_no;
int     n_inds_preference[number_of_cues];
CVec2D  centres[number_of_cues];

// output variables
int     cue_reached;
int     output_frequency;
double  path_length;
double  energy;
CVec2D  magnetisation;
CVec2D  centroid;
CVec2D* preference;
bool*   state;

std::string  filename;

// class vectors
spin*   agent;
cue*    CS;

// model functions
int main();
void RunGeneration();
void FlipSpins();
void CalculateSystemProperties(int spin_id);
void MoveAgents(int rep);
void CalculateEnergy(int spin_id);
void CalculateMagnetisation();
void SetupSimulation(double temp);
void SetupEnvironmentSymmetric();
void SetupEnvironmentAsymmetric();
void SetupEnvironmentFromCSV(const std::string& filename);
void SetupEnvironmentRandom();
void SetupSpins(double temp);
void ResetSetup(double x, double y);
void ResetStates();
CVec2D RandomBoundedPoint(double x, double y);
CVec2D RandomPolarPoint(double rmin, double rmax);
double GetProbability(double x, double mu, double sigma);

void GenerationOutput(int rep);

#endif /* parameteres_h */
