/*
 * map.h
 */

#ifndef MAP_H_
#define MAP_H_

class Map {
public:
	
	struct single_landmark_s{

		std::string id_i ; // Landmark ID
		double x_i; // Landmark x-position in the map (global coordinates)
        double y_i; // Landmark y-position in the map (global coordinates)
        double theta_i; //Landmark theta in the map (global coordinates)
	};

	std::vector<single_landmark_s> landmark_list ; // List of landmarks in the map

};



#endif /* MAP_H_ */
