#ifndef __MAKE_BUILDING_BLOCKS__
#define __MAKE_BUILDING_BLOCKS__
#include "buildingBlocks.hpp"
#include <vector>

BuildingBlocks makeBlocksCorridors(int nb_coll,
                                   std::vector<double>& corridor_1,
                                   std::vector<double>& corridor_2);

BuildingBlocks makeBlocksVirtualObstacles();

BuildingBlocks makeSimpleBlocks();

BuildingBlocks makeSimpleBlocksFreeTime();

BuildingBlocks makeBlocksManyObstacles();

BuildingBlocks makeBlocksMoonlander();

BuildingBlocks makeBlocksDiscontinuous(std::vector<double>& x0, 
                                       std::vector<double>& xf);

#endif