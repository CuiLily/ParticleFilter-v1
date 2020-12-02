/*
 * particle_filter.cpp
 */
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    num_particles = 729; //set to number of files in observation directory

    weights.resize(num_particles);
    particles.resize(num_particles);

    double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
    std_x = std[0];
    std_y = std[1];
    std_theta = std[2];

    // Normal distribution for x, y and theta
    normal_distribution<double> dist_x(x, std_x); // mean is centered around the new measurement
    normal_distribution<double> dist_y(y, std_y);
    normal_distribution<double> dist_theta(theta, std_theta);

    default_random_engine gen; //http://www.cplusplus.com/reference/random/default_random_engine/

    // create particles and set their values
    for(int i=0; i<num_particles; ++i){
        Particle p;
        p.id = i;
        p.x = dist_x(gen); // take a random value from the Gaussian Normal distribution and update the attribute
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
        p.weight = 1;

        particles[i] = p;
        weights[i] = p.weight;
    }
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
    std_x = std_pos[0];
    std_y = std_pos[1];
    std_theta = std_pos[2];

    default_random_engine gen;

    for(int i=0; i<num_particles; ++i){
        Particle *p = &particles[i]; // get address of particle to update

        // use the prediction equations from the Lesson 14
        double new_x = p->x + (velocity/yaw_rate) * (sin(p->theta + yaw_rate*delta_t) - sin(p->theta));
        double new_y = p->y + (velocity/yaw_rate) * (cos(p->theta) - cos(p->theta + yaw_rate*delta_t));
        double new_theta = p->theta + (yaw_rate*delta_t);

        // add Gaussian Noise to each measurement
        // Normal distribution for x, y and theta
        normal_distribution<double> dist_x(new_x, std_x);
        normal_distribution<double> dist_y(new_y, std_y);
        normal_distribution<double> dist_theta(new_theta, std_theta);

        // update the particle attributes
        p->x = dist_x(gen);
        p->y = dist_y(gen);
        p->theta = dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    for(auto pred : predicted){
      double dist_min = std::numeric_limits<double>::max();
      for(auto observation : observations){
        double distance = dist(observation.x, observation.y, pred.x, pred.y); // distance b/w obs and landmark
        if(distance < dist_min){
          observation.id = pred.id;
        }
        dist_min = distance;
      }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// Update the weights of each particle using a multi-variate Gaussian distribution.
	// NOTE: The observations are given in the image coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.

    double std_x = std_landmark[0];
    double std_y = std_landmark[1];
    double weights_sum = 0;

    for(int i=0; i<num_particles; ++i){
        Particle *p = &particles[i];
        double wt = 1.0;

        // convert observation from image to map's coordinate system
        for(int j=0; j<observations.size(); ++j){
            LandmarkObs current_obs = observations[j];
            //LandmarkObs transformed_obs;

            double fx = 8.6264757161231410e+02;
            double cx = 1.4419426569209147e+03;
            double fy = 8.6179474700008700e+02;
            double cy = 1.4373417306128667e+03;

            Eigen::Matrix3d K;
            K << fx, 0,  cx,
                  0, fy, cy,
                  0, 0,  1;

            Eigen::Vector3d px_homo(current_obs.x, current_obs.y, 1.0);
            Eigen::Vector3d f; //bearing vector

            double h = 1.77; //camera height: 1.77m
            Eigen::Vector3d c(p->x, p->y, h); //camera position

            // Rotation from camera frame to world
            Eigen::Matrix3d Rwc;
            Rwc << cos(p->theta), -sin(p->theta), 0.0,
                   sin(p->theta), cos(p->theta),  0.0,
                   0.0, 0.0, 1.0;

            // calculate bearing vector
            f = K.inverse() * px_homo;

            // depth: lambda
            double lambda = - h / f[2];
            // calculate s = (a, b, 0)
            double a = p->x + lambda * (f[0] * cos(p->theta) - f[1] * sin(p->theta));
            double b = p->y + lambda * (f[0] * sin(p->theta) + f[1] * cos(p->theta));

            Eigen::Vector3d s(a, b, 0.0); //center of parking lot code in the world

            // find the closest neighbor in the map list
            Map::single_landmark_s landmark;
            double distance_min = std::numeric_limits<double>::max();

            for(int k=0; k<map_landmarks.landmark_list.size(); ++k){
                Map::single_landmark_s cur_l = map_landmarks.landmark_list[k];
                double distance = dist(a, b, cur_l.x_i, cur_l.y_i);
                if(distance < distance_min){
                    distance_min = distance;
                    landmark = cur_l;
                }
            }


//            transformed_obs.x = (current_obs.x * cos(p->theta)) - (current_obs.y * sin(p->theta)) + p->x;
//            transformed_obs.y = (current_obs.x * sin(p->theta)) + (current_obs.y * cos(p->theta)) + p->y;
//            transformed_obs.id = current_obs.id;

//            // find the predicted measurement that is closest to each observed measurement and assign
//            // the observed measurement to this particular landmark
//            Map::single_landmark_s landmark;
//            double distance_min = std::numeric_limits<double>::max();
//
//            for(int k=0; k<map_landmarks.landmark_list.size(); ++k){
//                Map::single_landmark_s cur_l = map_landmarks.landmark_list[k];
//                double distance = dist(transformed_obs.x, transformed_obs.y, cur_l.x_i, cur_l.y_i);
//                if(distance < distance_min){
//                    distance_min = distance;
//                    landmark = cur_l;
//                }
//            }

            // update weights using Multivariate Gaussian Distribution
            // equation given in Transformations and Associations Quiz
            double num = exp(-0.5 * (pow((a - landmark.x_i), 2) / pow(std_x, 2) + pow((b - landmark.y_i), 2) / pow(std_y, 2)));
            double denom = 2 * M_PI * std_x * std_y;
            wt *= num/denom;
        }
        weights_sum += wt;
        p->weight = wt;
    }
    // normalize weights to bring them in (0, 1]
    for (int i = 0; i < num_particles; i++) {
        Particle *p = &particles[i];
        p->weight /= weights_sum;
        weights[i] = p->weight;
    }
}

void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    default_random_engine gen;

    // Random integers on the [0, n) range
    // the probability of each individual integer is its weight of the divided by the sum of all weights.
    discrete_distribution<int> distribution(weights.begin(), weights.end());
    vector<Particle> resampled_particles;

    for (int i = 0; i < num_particles; i++){
        resampled_particles.push_back(particles[distribution(gen)]);
    }

    particles = resampled_particles;

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
