/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

using namespace std;
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  
  num_particles = 100;  // TODO: Set the number of particles
  
  // Defination of sensor noise: Refer particle_filter.h
 
    // 1. declare std deviation
  double std_x = std [0];
  double std_y = std [1];
  double std_theta = std [2];
  
    // 2. Create Gaussian distribution 
  normal_distribution <double> x_init(x, std_x);
  normal_distribution <double> y_init(y, std_y);
  normal_distribution <double> theta_init(theta, std_theta);
  
    // 3. Create normal distribution with mean value as GPS values
  for (int i = 0; i < num_particles; ++i){
    
    Particle p;
    p.id = i;
    p.x = x_init(gen);
    p.y = y_init(gen);
    p.theta = theta_init(gen);
    p.weight = 1.0; // Initialize 1.0 weight to all particles
    
    particles.push_back(p);
  }
  
  is_initialized = true;
  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  // Define sensor noise and calculate new state
  
    // 1. Declare std deviation
  double std_x = std_pos [0];
  double std_y = std_pos [1];
  double std_theta = std_pos [2];  
  
    // 2. Create Gaussian distribution 
  normal_distribution <double> x_pred(0, std_x);
  normal_distribution <double> y_pred(0, std_y);
  normal_distribution <double> theta_pred(0, std_theta); 
  
    // 3. Calculate new state
  for (int i = 0; i < num_particles; ++i){
    double theta = particles[i].theta;
    
    if (fabs(yaw_rate) < 0.000001){ // Yaw rate is not changing or changing very minimally 
      particles[i].x += velocity * delta_t * cos(theta);
      particles[i].y += velocity * delta_t * sin(theta); 
    }
    else {
      particles[i].x += velocity / yaw_rate * (sin(theta + yaw_rate * delta_t) - sin(theta));
      particles[i].y += velocity / yaw_rate * (cos(theta) - cos(theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }
    
     // 4. Adding noise
    particles[i].x += x_pred(gen);
    particles[i].y += y_pred(gen);
    particles[i].theta += theta_pred(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  for (unsigned int i = 0; i < observations.size(); ++i){ // Each observation
    
    // 1. Current obs
    LandmarkObs obs = observations [i];
    
    // 2. Initialize min distance as big number
    double min_distance = numeric_limits<double>::max();
    
    // 3. Initialize map id a negative number
    int map_id = -1;
    
    for (unsigned int j = 0; j < predicted.size(); ++j){ // Each prediction
      LandmarkObs pred = predicted [j];
      
      // 4. Distnace between obs and predicted landmarks
      double distance = dist(obs.x, obs.y, pred.x, pred.y);
      
      if (distance < min_distance){
        min_distance = distance;
        map_id = pred.id;
      }
    }
    
    // 5. Set obs id to nearest id
    observations[i].id = map_id; 
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  // 1. Each particle
  for (int i = 0; i < num_particles; ++i){
    
    // 2. Get particle coordinates
    double p_x = particles [i].x;
    double p_y = particles [i].y;
    double p_theta = particles [i].theta;
    
    vector<LandmarkObs> predictions;
    
    // 3. Vector to hold predicted landmark location within sensor range
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j){
      float landmark_x = map_landmarks.landmark_list[j].x_f;
      float landmark_y = map_landmarks.landmark_list[j].y_f;
      int landmark_id = map_landmarks.landmark_list[j].id_i;
      
      // 4. Consider landmarks within the sensor range considering the rectangular region aroung particle (computationally faster than circular region consideration)
      if (fabs (landmark_x - p_x) <= sensor_range && fabs (landmark_y - p_y) <= sensor_range){
        
        // 5. Add this landmark to prediction
        predictions.push_back(LandmarkObs{landmark_id, landmark_x, landmark_y}); // Different 
      }
    }    
    // 5. Observations transformed from vehicle coordinates to map coordinates
    vector<LandmarkObs>transformed;
    for (unsigned int k = 0; k < observations.size(); ++k){
      double transform_x = cos(p_theta)*observations[k].x - sin(p_theta)*observations[k].y + p_x;
      double transform_y = sin(p_theta)*observations[k].x + cos(p_theta)*observations[k].y + p_y;
      transformed.push_back(LandmarkObs{observations[k].id, transform_x, transform_y});
    }
    
    // 6. Perform data association
    dataAssociation(predictions, transformed);
    
    particles[i].weight = 1.0;
    
    // 7. Placeholder for observations and predictions    associated with it 
    for (unsigned int l = 0; l < transformed.size(); ++l){
      double obs_x, obs_y, pred_x, pred_y;
      obs_x = transformed[l].x;
      obs_y = transformed[l].y;
      int associated_pred = transformed[l].id;
      
      for (unsigned int m = 0; m < predictions.size(); ++m){
        if (predictions[m].id == associated_pred){
          
          pred_x = predictions[m].x;
          pred_y = predictions[m].y;
        }          
      }
      
      // Calculate weights for the observation with multivariate Gaussian
      double s_x = std_landmark[0];
      double s_y = std_landmark[1];
      double obs_w = (1/(2*M_PI*s_x*s_y)) * exp(-(pow(pred_x-obs_x,2)/(2*pow(s_x, 2)) + (pow(pred_y-obs_y,2)/(2*pow(s_y, 2)))));
      
      particles[i].weight *= obs_w; 
    }
  }  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  vector<Particle> new_particles;
  
  // 1. Gather all current weights
  vector<double> weights;
  for (int i = 0; i < num_particles; ++i){
    weights.push_back(particles[i].weight);
  }
  
  // 2. Generate random index for resampling wheel
  uniform_int_distribution<int>int_dist(0, num_particles - 1);
  auto index = int_dist(gen);
  
  // 3. Get max weight
  double max_weight = *max_element(weights.begin(), weights.end());
  
  // 4. Random uniform distribution
  uniform_real_distribution<double>real_dist(0.0, max_weight);
  
  double beta = 0.0;
  
  // 5. Resample wheel spin
  for (int j = 0; j < num_particles; ++j){
    beta += real_dist(gen) * 2.0;
    while (beta > weights[index]){
      beta -= weights[index];
      index = (index + 1) % num_particles;      
    }
    new_particles.push_back(particles[index]);          
  }
  particles = new_particles;  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}