/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	num_particles = 10;
	default_random_engine gen;
  // Extracting standard deviations
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  // Creating normal distributions
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);


  for (int i = 0; i < num_particles; i++) {

    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    particles.push_back(particle);
	weights.push_back(particle.weight);
	}

  is_initialized = true;

  

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
  	double std_x = std_pos[0];
  	double std_y = std_pos[1];
  	double std_theta = std_pos[2];

  	

	fill(weights.begin(),weights.end(),1);
  for (int i = 0; i < num_particles; i++) {

		double theta = particles[i].theta;

		if ( fabs(yaw_rate) < 0.001 ) { 
		particles[i].x += velocity * delta_t * cos( theta );
		particles[i].y += velocity * delta_t * sin( theta );
		
		} else {
		particles[i].x += velocity / yaw_rate * ( sin( theta + yaw_rate * delta_t ) - sin( theta ) );
		particles[i].y += velocity / yaw_rate * ( cos( theta ) - cos( theta + yaw_rate * delta_t ) );
		particles[i].theta += yaw_rate * delta_t;
		}

		normal_distribution<double> dist_x(particles[i].x, std_x);
  		normal_distribution<double> dist_y(particles[i].y, std_y);
  		normal_distribution<double> dist_theta(particles[i].theta, std_theta);
		  
		// Adding noise.
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
		particles[i].weight = 1.0;

	}
	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	for(int i = 0; i < observations.size(); ++i) {
		LandmarkObs obs = observations[i];
		double min_dist = INFINITY;
		int closest_pid = -1;
		
		for(int j = 0; j < predicted.size(); ++j) {
		LandmarkObs pred = predicted[j];
		//double dist = dist(landmark_x,landmark_y,px,py);
		double distance = dist(obs.x, obs.y, pred.x, pred.y);

		if (distance < min_dist) {
			min_dist = distance;
			closest_pid = j;
			}
		}
		cout << min_dist <<" ";	
		observations[i].id = closest_pid;
  }
  

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	
	double sig_x = std_landmark[0] * std_landmark[0];
  	double sig_y = std_landmark[1] * std_landmark[1];
  	double normalizer = 2.0 * M_PI * std_landmark[0] * std_landmark[1];

	for(int i=0; i<num_particles; i++){
		double px = particles[i].x;
		double py = particles[i].y;
		double theta = particles[i].theta;
		
			
		double total_prob = 1.0;
	
		for(int k=0;k < observations.size();k++){

			double x_observ = observations[k].x;
			double y_observ = observations[k].y;
			int id_observ = observations[k].id;

			// transform to map coordinate
			double x_map= px + (cos(theta) * x_observ) - (sin(theta) * y_observ);
			double y_map= py + (sin(theta) * x_observ) + (cos(theta) * y_observ);

			double min_dist = 99999.99;
			int closest_pid = -1;
		
			for(int j = 0; j < map_landmarks.landmark_list.size(); ++j) {
				double distance = dist(x_map,y_map, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);

				if (distance < min_dist) {
					min_dist = distance;
					closest_pid = j;
				}
			}

			double dx = x_map - map_landmarks.landmark_list[closest_pid].x_f;
			double dy = y_map - map_landmarks.landmark_list[closest_pid].y_f;
			double power = (dx * dx) / (2 * sig_x) + (dy * dy) / (2 * sig_y);
			total_prob *= exp(-power) / normalizer;	
		}
		particles[i].weight = total_prob;
		weights[i] = total_prob;
		}
		  
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution


	vector<Particle> new_particles;

    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> distribution(weights.begin(), weights.end());

    for(int i = 0; i < num_particles; i++){
        Particle p = particles[distribution(gen)];
        new_particles.push_back(p);
    }
    particles = new_particles;

}



Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
