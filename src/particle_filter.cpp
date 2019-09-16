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
using std::normal_distribution;

  // Set the random number generator
static  std::default_random_engine r_generator; 

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	/**
	* TODO: Set the number of particles. Initialize all particles to 
	*   first position (based on estimates of x, y, theta and their uncertainties
	*   from GPS) and all weights to 1. 
	* TODO: Add random Gaussian noise to each particle.
	* NOTE: Consult particle_filter.h for more information about this method 
	*   (and others in this file).
	*/
	if(is_initialized==false) //The initialisation should run only once.
	{ 
		num_particles = 2000;  // TODO: Set the number of particles
		//Create the normal distribution for x, y and theta
		normal_distribution<double> dist_x(x, std[0]);
		normal_distribution<double> dist_y(y, std[1]);
		normal_distribution<double> dist_theta(theta, std[2]);
		// Generate the particles
		for(int i=0;i<num_particles;i++)
		{
			Particle p; //Create a new Particle object and fill the values
			p.id = i;
			p.x = dist_x(r_generator);
			p.y = dist_y(r_generator);
			p.theta = dist_theta(r_generator);
			p.weight = 1.0;
			
			//Add the p to the end of the 'particles' vector
			particles.push_back(p); 
		}
	// Change the status of init indicator
	is_initialized = true;
	}
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

	for(int i=0;i<num_particles;i++)
	{
		Particle& p = particles[i];
		if(fabs(yaw_rate) > 0.00001)
		{
			p.x += (velocity/yaw_rate)*(sin(p.theta+yaw_rate*delta_t)-sin(p.theta));
			p.y += (velocity/yaw_rate)*(cos(p.theta)-cos(p.theta+yaw_rate*delta_t));
			p.theta += yaw_rate*delta_t;
		}
		else
		{
			p.x += velocity*delta_t*cos(p.theta);
			p.y += velocity*delta_t*sin(p.theta);
			// The theta is unchanged			
		}
		
		//Finally add noise
		normal_distribution<double> dist_x(p.x, std_pos[0]);
		normal_distribution<double> dist_y(p.y, std_pos[1]);
		normal_distribution<double> dist_theta(p.theta, std_pos[2]);

		p.x = dist_x(r_generator);
		p.y = dist_y(r_generator);
		p.theta = dist_theta(r_generator);	
		
		//std::cout<<"New prediction:" << particles[4].x;
		
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
   
   // Let's go through the observation points
   for(int i=0;i<observations.size();i++)
   {
	   // We need a big number for the initial distance
	   double min_distance = 999999999;
	   
	   // I set the detected landmark id to 0
	   // It will be overwritten.
	   int map_lm_id = 0;
	   
	   //Lets go through the predicted points
	   for(int j=0; j<predicted.size(); j++)
	   {
			// Calculate the current distance square, square root needs more time.
			double curr_distance = dist_sq(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
			
			if(curr_distance<min_distance)
			{
				map_lm_id = predicted[j].id;
				min_distance = curr_distance;
			}
	   }	   
	   observations[i].id = map_lm_id;
   }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   vector<LandmarkObs> &observations, 
                                   Map &map_landmarks) {
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
  
   // This shouldn't be reperead
      double sensor_range_sq = sensor_range*sensor_range;
	  double lm_x, lm_y;
   
   // Go through all particles
   for(int i=0; i< num_particles; i++)
   {
		// Generate a list for landmarks within sensor range
		vector<LandmarkObs> curr_landmarks;
		Particle& p = particles[i];
	   
		for(int j=0; j<map_landmarks.landmark_list.size(); j++)
		{
		   Map::single_landmark_s& l = map_landmarks.landmark_list[j];
		   if(dist_sq(p.x,p.y,l.x_f,l.y_f) <= sensor_range_sq)
		   {
			   curr_landmarks.push_back(LandmarkObs{ l.id_i, l.x_f, l.y_f });
			}	
		}
		
		// Transform the observation coordinates to map coordinate system
		vector<LandmarkObs> tr_observations;
		
		for(int j=0; j<observations.size(); j++)
		{
			LandmarkObs& o = observations[j];
			double map_x = p.x + cos(p.theta)*o.x - sin(p.theta)*o.y;
			double map_y = p.y + sin(p.theta)*o.x + cos(p.theta)*o.y;
			
			// Load the calculated result to the tr_observations vector.
			tr_observations.push_back(LandmarkObs{o.id, map_x, map_y});
		}
				
		// Associate the sensed landmarks(within sensor range) and the transformed observations
		dataAssociation(curr_landmarks,tr_observations);
		
		// Reset weight 
		p.weight = 1.0;
		
		for(int k=0; k<tr_observations.size(); k++)
		{
			LandmarkObs& t = tr_observations[k];
			int curr_lm_id = t.id;

			for(int l=0; l<curr_landmarks.size(); l++)
			{
								if(curr_lm_id == curr_landmarks[l].id)
				{
					lm_x = curr_landmarks[l].x;
					lm_y = curr_landmarks[l].y;
					break;					
				}
			}

			//Weight for this observation with multivariate Gaussian
			double obs_weight = ( 1/(2*M_PI*std_landmark[0]*std_landmark[1])) 
							* exp( -( pow(t.x-lm_x,2)/(2*pow(std_landmark[0], 2)) 
							+ (pow(t.y-lm_y,2)/(2*pow(std_landmark[1], 2))) ) );

			//Product of this obersvation weight with total observations weight
			particles[i].weight *= obs_weight;	
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


	double beta_max=0;

	vector<double> weights;
	double mw = std::numeric_limits<double>::min();
	for(int i=0; i<num_particles;i++)
	{
	   double curr_weight = particles[i].weight;
	   weights.push_back(curr_weight);
	   beta_max += curr_weight;
	   if(curr_weight>mw){
		   mw=curr_weight;
	   }
	}
  
	// Resample
	/**
	I implemented two resampling technic:
	1.) Residual
	2.) Systematic
	The second one gave tha same accuracy, but it's slightly quicker.
	Ref: https://pdfs.semanticscholar.org/8fcb/1ca95d38a2e014470856c1dfa21212aa90a5.pdf
	*/
	
	//std::uniform_real_distribution<double> distDouble(0.0, mw); //For residual
	std::uniform_real_distribution<double> distDouble(0.0, beta_max); //For systematic
	std::uniform_int_distribution<int> distInt(0, num_particles - 1);
	int index = distInt(r_generator)% num_particles;
   
	vector<Particle> new_particles;
	double beta = distDouble(r_generator);
	double step = beta_max/num_particles;
	
	for(int i=0; i<particles.size(); i++)
	{
		//beta += distDouble(r_generator) * 2.0; //For residual
		beta += step; //For systematic (slightly quicker)
		while(beta > weights[index]) 
		{
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