//
//  parameteres.h
//  multi-choice_decision_geometry
//
//  Created by Vivek Hari Sridhar on 18/09/17.
//  Copyright © 2017 Vivek Hari Sridhar. All rights reserved.
//

#ifndef parameteres_h
#define parameteres_h

#include "spin.h"
#include <fstream>

int     timestep_number;    // timestep number
int     arena_size;
int     total_agents;
int     cue_reached;
int     trial_time;
int     reset_no;
int     grid_cols;          // create a lattice to ensure neighbour consistency
int     grid_rows;
int     interaction_radius;

double  spin_energy;
CVec2D  system_magnetisation;

double  start_dist;
double  overall_angle;
double  dist_thresh;
double  max_angle;
double  path_length;
double  avg_speed;

bool    rep_done;
bool    symmetric;

CVec2D  arena_centre;
CVec2D  centroid;
CVec2D  centres[number_of_cues];
int     n_inds_preference[number_of_cues];

spin*   agent;
cue*    CS;

int main();
void RunGeneration();
void FlipSpins();
void CalculateSpinProperties(int spin_id);
void MoveAgents();
void SetupSimulation(double temp);
void SetupEnvironmentSymmetric();
void SetupEnvironmentAsymmetric();
void SetupSpins(double temp);
void ResetSetup();
CVec2D RandomBoundedPoint();
void GenerationalOutput(double temp, int rep);

void Graphics();
void GraphicsWriter(cv::VideoWriter& video_writer, int& timestep_number, const int& num_timesteps);

#endif /* parameteres_h */
