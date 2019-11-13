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

double  system_energy;
CVec2D  system_magnetisation;

double  left_right_dist;
double  start_dist;
double  dist_thresh;
double  max_angle;
double  hat_width;

double  path_length;

bool    rep_done;
bool    symmetric;
bool    distance;

CVec2D  arena_centre;
CVec2D  centroid;
CVec2D  centres[number_of_cues];
int     n_inds_preference[number_of_cues];

spin*   agent;
cue*    CS;

int main();
void RunGeneration();
void FlipSpins();
void CalculateSystemProperties(int spin_id);
void CalculateSpinProperties(double& energy, int spin_id);
void MoveAgents(int rep);
void SetupSimulation(double temp);
void SetupEnvironmentSymmetric();
void SetupEnvironmentAsymmetric();
void SetupEnvironmentDistances();
void SetupSpins(double temp);
void ResetSetup(double x, double y);
CVec2D RandomBoundedPoint(double x, double y);
void GenerationalOutput(int rep);

void Graphics();

#endif /* parameteres_h */
