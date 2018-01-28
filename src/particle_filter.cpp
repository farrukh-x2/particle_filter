/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Framework: Tiffany Huang
 *      Completed: Farrukh Ali
 *
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
#include <tuple>

#include "particle_filter.h"

using namespace std;

//Functions and variables not defined in the .h file, submissions only accepted with particle_filter.cpp file
std::default_random_engine gen;

void dataAssociationByMap(std::vector<LandmarkObs>& transformed_obs, Map map_landmarks);
std::tuple<double, double> transformedObs(double particle_x,
       double particle_y, double particle_theta, double obs_x, double obs_y);


void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
    num_particles = 25;
    is_initialized = true;

    // This creates a normal (Gaussian) distribution for x, y and theta.
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    for (int i = 0; i < num_particles; ++i) {

        Particle single_particle;
        //Sample these from normal distributions
        single_particle.x = dist_x(gen);
        single_particle.y = dist_y(gen);
        single_particle.theta = dist_theta(gen);
        single_particle.weight = 1.0;
        single_particle.id = i;

        particles.push_back(single_particle);
        weights.push_back(1.0);

    }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // Add measurements to each particle and add random Gaussian noise.

    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);

    if (fabs(yaw_rate) < 0.0001){
        yaw_rate = 0.0001;
    }

    for (int i = 0; i < num_particles; ++i) {

        double x0 = particles[i].x;
        double y0 = particles[i].y;
        double theta_0 = particles[i].theta;
        //velocity = 110;
        //yaw_rate = 3.142/8;
        double xf, yf, theta_f;

        xf = x0 + (velocity/yaw_rate)*(sin(theta_0 + (yaw_rate*delta_t))  - sin(theta_0) );
        yf = y0 + (velocity/yaw_rate)*(cos(theta_0) - cos(theta_0  + (yaw_rate*delta_t)) );
        theta_f = theta_0 + (yaw_rate*delta_t) ;

        particles[i].x = xf + dist_x(gen);
        particles[i].y = yf + dist_y(gen);
        particles[i].theta = theta_f + dist_theta(gen);

    }

}

void dataAssociationByMap(std::vector<LandmarkObs>& transformed_obs, Map map_landmarks) {
    //  Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
    double temp_dis;

    for (int i=0; i< transformed_obs.size(); ++i){
        double dis_landmark_obs = 999999;

        for(int j = 0; j < map_landmarks.landmark_list.size(); ++j)
        {
            temp_dis = dist(map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f,transformed_obs[i].x,transformed_obs[i].y);
            //cout<< "temp_dis " << temp_dis<<endl;
            if (temp_dis < dis_landmark_obs)
            {
                dis_landmark_obs = temp_dis;
                transformed_obs[i].id = map_landmarks.landmark_list[j].id_i;
            }
        }

    }

}

std::tuple<double, double> transformedObs(double particle_x,
       double particle_y, double particle_theta, double obs_x, double obs_y){
       //The observations are given in the VEHICLE'S coordinate system. The particles are located
       //according to the MAP'S coordinate system. Need to transform between the two systems.
       //This transformation requires both rotation AND translation (but no scaling).

    double transf_x, transf_y;
    transf_x = obs_x*cos(particle_theta) - obs_y*sin(particle_theta) + particle_x;
    transf_y = obs_x*sin(particle_theta) + obs_y*cos(particle_theta) + particle_y;

    return std::make_tuple(transf_x, transf_y);

}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
    //  Update the weights of each particle using a multi-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution

    double weights_sum = 0.0;
    for (int k = 0; k < particles.size(); ++k){

        vector<LandmarkObs> transformed_obs_vector;
        for(int i = 0; i < observations.size(); ++i)
        {
            LandmarkObs trans_obs;
            tie(trans_obs.x,trans_obs.y) = transformedObs(particles[k].x,particles[k].y,particles[k].theta,observations[i].x,observations[i].y);
            transformed_obs_vector.push_back(trans_obs);
        }

        dataAssociationByMap(transformed_obs_vector, map_landmarks );

        double t_prob = 1.0;
        for (int jj=0; jj< transformed_obs_vector.size(); ++jj){
        double x = transformed_obs_vector[jj].x; //0;
        double y = transformed_obs_vector[jj].y; //5;

        double ux = map_landmarks.landmark_list[transformed_obs_vector[jj].id-1].x_f;
        double uy = map_landmarks.landmark_list[transformed_obs_vector[jj].id-1].y_f;

        double A = ( x - ux )*( x - ux) / (2 * std_landmark[0]* std_landmark[0]);
        double B = ( y - uy )*( y - uy) / (2 * std_landmark[1]* std_landmark[1]);
        t_prob = t_prob * (1 / (2*M_PI*std_landmark[0]*std_landmark[1])) * exp(-1 * (A + B));

        }

        particles[k].weight = t_prob;
        weights[k] = t_prob;
        weights_sum += t_prob;
    }

    //Normalize the wieghts vector
    for(int i = 0; i < weights.size(); ++i){
                weights[i] /= weights_sum;
                //if(weights[i] > 0)
                 //   cout<<i<<":"<<weights[i]<<endl;
       }

}
void ParticleFilter::resample() {
    // Resample particles with replacement with probability proportional to their weight.

    std::discrete_distribution<int> dist(weights.begin(), weights.end());
    std::vector<Particle> rsmpld_prtcls(particles.size());
    for(int i = 0; i < num_particles; i++){
        int ind = dist(gen);
        rsmpld_prtcls[i] = particles[ind];
    }
    particles = rsmpld_prtcls;

}
//Functions for diagnostics and optional use
Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
