/*
 * helper_functions.h
 * Some helper functions for the 2D particle filter.
 */

#ifndef HELPER_FUNCTIONS_H_
#define HELPER_FUNCTIONS_H_

#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include "map.h"

using namespace std;
/*
 * Struct representing one position/control measurement.
 */
struct control_s {
	
	double velocity;	// Velocity [m/s]
	double yawrate;		// Yaw rate [rad/s]
};

/*
 * Struct representing one ground truth position.
 */
struct ground_truth {
	
	double x;		// Global vehicle x position [m]
	double y;		// Global vehicle y position
	double theta;	// Global vehicle yaw [rad]
};

/*
 * Struct representing one landmark observation measurement.
 */
struct LandmarkObs {

    int id;
    // position of center
    double x;			// Local (vehicle coordinates) x position of landmark observation [m]
    double y;			// Local (vehicle coordinates) y position of landmark observation [m]

	// probability for each candidates of each digit
	double d10, d11, d12, d13, d14, d15, d16, d17, d18 ,d19; // digit1
	double d20, d21, d22, d23, d24, d25, d26, d27, d28, d29; // digit2
	double d30, d31, d32, d33, d34, d35, d36, d37, d38, d39; // digit3
};

/*
 * Computes the Euclidean distance between two 2D points.
 * @param (x1,y1) x and y coordinates of first point
 * @param (x2,y2) x and y coordinates of second point
 * @output Euclidean distance between two 2D points
 */
inline double dist(double x1, double y1, double x2, double y2) {
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

inline double * getError(double gt_x, double gt_y, double gt_theta, double pf_x, double pf_y, double pf_theta) {
	static double error[3];
	error[0] = fabs(pf_x - gt_x);
	error[1] = fabs(pf_y - gt_y);
	error[2] = fabs(pf_theta - gt_theta);
	error[2] = fmod(error[2], 2.0 * M_PI);
	if (error[2] > M_PI) {
		error[2] = 2.0 * M_PI - error[2];
	}
	return error;
}

/* Reads map data from a file.
 * @param filename Name of file containing map data.
 * @output True if opening and reading file was successful
 */
inline bool read_map_data(std::string filename, Map& map) {

	// Get file of map:
	std::ifstream in_file_map(filename.c_str(),std::ifstream::in);//ifstream:输入文件流，用于从文件读取信息。
	// Return if we can't open the file.
	if (!in_file_map) {
		return false;
	}
    //ignore the 1st line of txt file
    std::string dummyLine;
    getline(in_file_map, dummyLine);
	
	// Declare single line of map file:
	std::string line_map;

	// Run over each single line:
	while(getline(in_file_map, line_map)){

		std::istringstream iss_map(line_map);

		// Declare landmark values and ID:
        std::string id_i;
        double landmark_x_i, landmark_y_i, landmark_theta_i;

		// Read data from current line to values::
        iss_map >> id_i;
		iss_map >> landmark_x_i; //把iss_map中的，以空格隔开的字符提取出来
		iss_map >> landmark_y_i;
		iss_map >> landmark_theta_i;

		cout<<"id_i: "<<id_i<<endl;
        cout << "landmark_x_i: " << landmark_x_i << endl;
        cout << "landmark_y_i: " << landmark_y_i << endl;
        cout << "landmark_theta_i: " << landmark_theta_i << endl;

		// Declare single_landmark: 是Map类的一个结构体
		Map::single_landmark_s single_landmark_temp; 

		// Set values 设置map的对象
		single_landmark_temp.id_i = id_i;
		single_landmark_temp.x_i  = landmark_x_i;
		single_landmark_temp.y_i  = landmark_y_i;
		single_landmark_temp.theta_i = landmark_theta_i;

		// Add to landmark list of map:添加到map中的landmark_list
		map.landmark_list.push_back(single_landmark_temp);
	}
	return true;
}

/* Reads control data from a file.
 * @param filename Name of file containing control measurements.
 * @output True if opening and reading file was successful
 */
inline bool read_control_data(std::string filename, std::vector<control_s>& position_meas) {

	// Get file of position measurements:
	std::ifstream in_file_pos(filename.c_str(),std::ifstream::in);
	// Return if we can't open the file.
	if (!in_file_pos) {
		return false;
	}
    //ignore the 1st line of txt file
    std::string dummyLine;
    getline(in_file_pos, dummyLine);
	// Declare single line of position measurement file:
	std::string line_pos;

	// Run over each single line:
	while(getline(in_file_pos, line_pos)){

		std::istringstream iss_pos(line_pos);

		// Declare position values:
		double velocity, yawrate;

		// Declare single control measurement: control_s是结构体
		control_s meas; 

		//read data from line to values:

		iss_pos >> velocity; //把iss_pos里的字符提取出来
		iss_pos >> yawrate;

        cout << "velocity: " << velocity << endl;
        cout << "yaw rate: " << yawrate << endl;

		// Set values
		meas.velocity = velocity; //赋值给结构体的这个对象
		meas.yawrate = yawrate;

		// Add to list of control measurements:
		position_meas.push_back(meas); //把这个测量值放进position_meas这个vector
	}
	return true;
}

/* Reads ground truth data from a file.
 * @param filename Name of file containing ground truth.
 * @output True if opening and reading file was successful
 */
inline bool read_gt_data(std::string filename, std::vector<ground_truth>& gt) {

	// Get file of position measurements:
	std::ifstream in_file_pos(filename.c_str(),std::ifstream::in);
	// Return if we can't open the file.
	if (!in_file_pos) {
		return false;
	}
    //ignore the 1st line of txt file
    std::string dummyLine;
    getline(in_file_pos, dummyLine);
	// Declare single line of position measurement file:
	std::string line_pos;

	// Run over each single line:
	while(getline(in_file_pos, line_pos)){

		std::istringstream iss_pos(line_pos);

		// Declare position values:
		double x, y, azimuth;//方位角（就是围绕z轴转的角度）

		// Declare single ground truth:
		ground_truth single_gt; //一个ground truth

		//read data from line to values:
		iss_pos >> x;
		iss_pos >> y;
		iss_pos >> azimuth;

        cout << "x: " << x << endl;
        cout << "y: " << y << endl;
        cout<<"azimuth: " << azimuth <<endl;

		// Set values
		single_gt.x = x;
		single_gt.y = y;
		single_gt.theta = azimuth;

		// Add to list of control measurements and ground truth:
		gt.push_back(single_gt);
	}
	return true;
}

/* Reads landmark observation data from a file.
 * @param filename Name of file containing landmark observation measurements.
 * @output True if opening and reading file was successful
 */
inline bool read_landmark_data(std::string filename, std::vector<LandmarkObs>& observations) {

	// Get file of landmark measurements:
	std::ifstream in_file_obs(filename.c_str(),std::ifstream::in);
	// Return if we can't open the file.
	if (!in_file_obs) {
		return false;
	}
    //ignore the 1st line of txt file
    std::string dummyLine;
    getline(in_file_obs, dummyLine);

	// Declare single line of landmark measurement file:
	std::string line_obs;

	// Run over each single line:
	while(getline(in_file_obs, line_obs)){

		std::istringstream iss_obs(line_obs);

		// Declare position values:
		double local_x, local_y;
        double d10, d11, d12, d13, d14, d15, d16, d17, d18 ,d19; // digit1
        double d20, d21, d22, d23, d24, d25, d26, d27, d28, d29; // digit2
        double d30, d31, d32, d33, d34, d35, d36, d37, d38, d39; // digit3

		//read data from line to values:
		iss_obs >> local_x;
		iss_obs >> local_y;
		// for digit1
		iss_obs >> d10;
        iss_obs >> d11;
        iss_obs >> d12;
        iss_obs >> d13;
        iss_obs >> d14;
        iss_obs >> d15;
        iss_obs >> d16;
        iss_obs >> d17;
        iss_obs >> d18;
        iss_obs >> d19;
        // for digit2
        iss_obs >> d20;
        iss_obs >> d21;
        iss_obs >> d22;
        iss_obs >> d23;
        iss_obs >> d24;
        iss_obs >> d25;
        iss_obs >> d26;
        iss_obs >> d27;
        iss_obs >> d28;
        iss_obs >> d29;
        // for digit
        iss_obs >> d30;
        iss_obs >> d31;
        iss_obs >> d32;
        iss_obs >> d33;
        iss_obs >> d34;
        iss_obs >> d35;
        iss_obs >> d36;
        iss_obs >> d37;
        iss_obs >> d38;
        iss_obs >> d39;

        cout << "local_x: " << local_x << endl;
        cout << "local_y: " << local_y << endl;

		// Declare single landmark measurement:
		LandmarkObs meas; //LandmarkObs是一个结构体

		// Set values
		meas.x = local_x;
		meas.y = local_y;
		//for digit2
		meas.d10 = d10;
        meas.d11 = d11;
        meas.d12 = d12;
        meas.d13 = d13;
        meas.d14 = d14;
        meas.d15 = d15;
        meas.d16 = d16;
        meas.d17 = d17;
        meas.d18 = d18;
        meas.d19 = d19;
        //for digit1
        meas.d20 = d20;
        meas.d21 = d21;
        meas.d22 = d22;
        meas.d23 = d23;
        meas.d24 = d24;
        meas.d25 = d25;
        meas.d26 = d26;
        meas.d27 = d27;
        meas.d28 = d28;
        meas.d29 = d29;
        //for digit3
        meas.d30 = d30;
        meas.d31 = d31;
        meas.d32 = d32;
        meas.d33 = d33;
        meas.d34 = d34;
        meas.d35 = d35;
        meas.d36 = d36;
        meas.d37 = d37;
        meas.d38 = d38;
        meas.d39 = d39;

		// Add to list of control measurements:
		observations.push_back(meas); //把meas存到observation这个vector里
	}
	return true;
}

#endif /* HELPER_FUNCTIONS_H_ */
