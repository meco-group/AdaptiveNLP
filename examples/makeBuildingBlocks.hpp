#ifndef __MAKE_BUILDING_BLOCKS__
#define __MAKE_BUILDING_BLOCKS__
#include "core.hpp"
#include <vector>

BuildingBlocks minimalBlocks();

BuildingBlocks makeBlocksCorridors(int nb_coll,
                                   std::vector<double>& corridor_1,
                                   std::vector<double>& corridor_2);

BuildingBlocks makeBlocksManyObstacles();

BuildingBlocks makeBlocksMoonlander();

#endif