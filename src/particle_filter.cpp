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

using std::normal_distribution;
using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[])
{
  num_particles = 50; // TODO: Set the number of particles
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::default_random_engine gen;
  double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
  // Set standard deviations for x, y, and theta
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];

  // This line creates a normal (Gaussian) distribution
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  // Print your initial values to the terminal.
  std::cout << "Initial value " << x << " " << y << " "
            << theta << std::endl;
  for (int j = 0; j < num_particles; ++j)
  {
    double sample_x, sample_y, sample_theta;
    Particle a_particle;
    // Sample from these normal distributions like this:
    //   sample_x = dist_x(gen);
    //   where "gen" is the random engine initialized earlier.
    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);
    a_particle.x = sample_x;
    a_particle.y = sample_y;
    a_particle.theta = sample_theta;
    a_particle.id = j;
    a_particle.weight = 1.0;
    // Print your samples to the terminal.
    std::cout << "Sample " << j + 1 << " " << a_particle.x << " " << a_particle.y << " "
              << a_particle.theta << std::endl;
    particles.push_back(a_particle);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate)
{
  /**
   * Add measurements to each particle and add random Gaussian noise.
   */
  std::default_random_engine gen;
  double x_f, y_f, theta_f;
  for (int j = 0; j < num_particles; ++j)
  {
    double prev_theta = particles.at(j).theta;
    theta_f = particles.at(j).theta + yaw_rate * delta_t;
    if (yaw_rate != 0.0f)
    {
      x_f = particles.at(j).x + velocity / yaw_rate * (sin(theta_f) - sin(prev_theta));
      y_f = particles.at(j).y + velocity / yaw_rate * (cos(prev_theta) - cos(theta_f));
    }
    else
    {
      x_f = particles.at(j).x + velocity * delta_t * cos(theta_f);
      y_f = particles.at(j).y + velocity * delta_t * sin(theta_f);
    }
    normal_distribution<double> dist_x(x_f, std_pos[0]);
    normal_distribution<double> dist_y(y_f, std_pos[1]);
    normal_distribution<double> dist_theta(theta_f, std_pos[2]);
    particles.at(j).x = dist_x(gen);
    particles.at(j).y = dist_y(gen);
    particles.at(j).theta = dist_theta(gen);
  }
}

vector<LandmarkObs> ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs> &observations)
{
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   */
  vector<double> distance_to_landmark(observations.size(), 10000);
  vector<LandmarkObs> seen_landmarks(observations.size());
  for (unsigned int i = 0; i < predicted.size(); ++i)
  {
    for (unsigned int j = 0; j < observations.size(); ++j)
    {
      double current_dist = sqrt(pow(observations.at(j).x - predicted.at(i).x, 2.0f) + pow(observations.at(j).y - predicted.at(i).y, 2.0f));
      if ((current_dist < distance_to_landmark.at(j)))
      {
        distance_to_landmark.at(j) = current_dist;
        observations.at(j).id = predicted.at(i).id;
        seen_landmarks.at(j) = predicted.at(i);
      }
    }
  }
  //for (unsigned int j = 0; j < observations.size(); ++j){
  //  std::cout << j << ")  obsX: "<< observations.at(j).x << " landX:  "<< seen_landmarks.at(j).x << "obsY: "<< observations.at(j).y << " landY:  "<< seen_landmarks.at(j).y  << std::endl;
  //}
  return seen_landmarks;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks)
{
  /**
   * Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   * So they must be transformed
   */

  //Landmarks vector
  vector<LandmarkObs> landmarks;
  for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j){
    LandmarkObs landmark;
    landmark.id = map_landmarks.landmark_list.at(j).id_i;
    landmark.x = map_landmarks.landmark_list.at(j).x_f;
    landmark.y = map_landmarks.landmark_list.at(j).y_f;
    landmarks.push_back(landmark);
  }
  weights.clear();
  vector<LandmarkObs> map_observations;
  for(unsigned int i = 0; i < particles.size(); ++i){ 
    //std::cout << "Particle "<< i << std::endl;
    //Transform observations to map global frame
    LandmarkObs obs_map;
    map_observations.clear();
    for (unsigned int j = 0; j < observations.size(); ++j){
      obs_map.x = observations.at(j).x * cos(particles.at(i).theta) - observations.at(j).y * sin(particles.at(i).theta) + particles.at(i).x;
      obs_map.y = observations.at(j).x * sin(particles.at(i).theta) + observations.at(j).y * cos(particles.at(i).theta) + particles.at(i).y;
      obs_map.id = 0;
      double distance_to_landmark = sqrt(observations.at(j).x*observations.at(j).x + observations.at(j).y*observations.at(j).y);
      if(distance_to_landmark < sensor_range)
        map_observations.push_back(obs_map);
    }
    
    //Associate observations to measurements
    vector<LandmarkObs> seen_landmarks = dataAssociation(landmarks, map_observations);
    // Calculate probability for each measurement using multivariate gaussian distribution
    particles.at(i).weight = 1.0f;
    for (unsigned int j = 0; j < map_observations.size(); ++j){
      double p1 = pow(map_observations.at(j).x - seen_landmarks.at(j).x, 2.0f) / (2*pow(std_landmark[0],2.0f));
      double p2 = pow(map_observations.at(j).y - seen_landmarks.at(j).y, 2.0f) / (2*pow(std_landmark[1],2.0f));
      double p_x_y = exp(-(p1 + p2)) / (2*std_landmark[0]*std_landmark[1]);
      
      particles.at(i).weight = particles.at(i).weight*p_x_y;
      //std::cout << "Observation "<< j << " weight " << p_x_y << std::endl;
    }
    weights.push_back(particles.at(i).weight);
    //std::cout << "Particle "<< i << " weight "<< particles.at(i).weight << std::endl;
  }
}

void ParticleFilter::resample()
{
  /**
   * Resample particles with replacement with probability proportional 
   *   to their weight. 
   */
  std::default_random_engine gen;
  std::discrete_distribution<> d_distribution(weights.begin(),weights.end());
  std::vector<Particle> new_particles;
  
  for (int i = 0; i < num_particles; i++){
    int index = d_distribution(gen);
    new_particles.push_back(particles[index]);
  }
  particles = new_particles;

}

void ParticleFilter::SetAssociations(Particle &particle,
                                     const vector<int> &associations,
                                     const vector<double> &sense_x,
                                     const vector<double> &sense_y)
{
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord)
{
  vector<double> v;

  if (coord == "X")
  {
    v = best.sense_x;
  }
  else
  {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}