#include <explore/frontier_search.h>

#include <mutex>

#include <costmap_2d/cost_values.h>
#include <costmap_2d/costmap_2d.h>
#include <geometry_msgs/Point.h>
#include <explore/costmap_tools.h>

#define PI 3.14159265

namespace frontier_exploration
{
using costmap_2d::LETHAL_OBSTACLE;
using costmap_2d::NO_INFORMATION;
using costmap_2d::FREE_SPACE;

FrontierSearch::FrontierSearch(costmap_2d::Costmap2D* costmap,
                               double potential_scale, double gain_scale,
                               double min_frontier_size)
  : costmap_(costmap)
  , potential_scale_(potential_scale)
  , gain_scale_(gain_scale)
  , min_frontier_size_(min_frontier_size)
{
}

std::vector<Frontier> FrontierSearch::searchFrom(geometry_msgs::Point position)
{
  std::vector<Frontier> frontier_list;

  // Sanity check that robot is inside costmap bounds before searching
  unsigned int mx, my;
  if (!costmap_->worldToMap(position.x, position.y, mx, my)) {
    ROS_ERROR("Robot out of costmap bounds, cannot search for frontiers");
    return frontier_list;
  }

  // make sure map is consistent and locked for duration of search
  std::lock_guard<costmap_2d::Costmap2D::mutex_t> lock(*(costmap_->getMutex()));

  map_ = costmap_->getCharMap();
  size_x_ = costmap_->getSizeInCellsX();
  size_y_ = costmap_->getSizeInCellsY();

  // initialize flag arrays to keep track of visited and frontier cells
  std::vector<bool> frontier_flag(size_x_ * size_y_, false);
  std::vector<bool> visited_flag(size_x_ * size_y_, false);

  // initialize breadth first search
  std::queue<unsigned int> bfs;

  // find closest clear cell to start search
  unsigned int clear, pos = costmap_->getIndex(mx, my);
  if (nearestCell(clear, pos, FREE_SPACE, *costmap_)) {
    bfs.push(clear);
  } else {
    bfs.push(pos);
    ROS_WARN("Could not find nearby clear cell to start search");
  }
  visited_flag[bfs.front()] = true;

  while (!bfs.empty()) {
    unsigned int idx = bfs.front();
    bfs.pop();

    // iterate over 4-connected neighbourhood
    for (unsigned nbr : nhood4(idx, *costmap_)) {
      // add to queue all free, unvisited cells, use descending search in case
      // initialized on non-free cell
      if (map_[nbr] <= map_[idx] && !visited_flag[nbr]) {
        visited_flag[nbr] = true;
        bfs.push(nbr);
        // check if cell is new frontier cell (unvisited, NO_INFORMATION, free
        // neighbour)
      } else if (isNewFrontierCell(nbr, frontier_flag)) {
        frontier_flag[nbr] = true;
        Frontier new_frontier = buildNewFrontier(nbr, pos, frontier_flag);
        if (new_frontier.size * costmap_->getResolution() >=
            min_frontier_size_) {
          frontier_list.push_back(new_frontier);
        }
      }
    }
  }

  // set costs of frontiers
  for (auto& frontier : frontier_list) {
    frontier.cost = frontierCost(frontier);
  }
  std::sort(
      frontier_list.begin(), frontier_list.end(),
      [](const Frontier& f1, const Frontier& f2) { return f1.cost < f2.cost; });

  return frontier_list;
}

Frontier FrontierSearch::buildNewFrontier(unsigned int initial_cell,
                                          unsigned int reference,
                                          std::vector<bool>& frontier_flag)
{
  // initialize frontier structure
  Frontier output;
  output.centroid.x = 0;
  output.centroid.y = 0;
  output.size = 1;
  output.min_distance = std::numeric_limits<double>::infinity();
  // [chad] add orientation info to ensure the racecar facing toward frontiers
  // theta is the orientation facing toward unknown area
  output.theta = 0.0;


  // record initial contact point for frontier
  unsigned int ix, iy;
  costmap_->indexToCells(initial_cell, ix, iy);
  costmap_->mapToWorld(ix, iy, output.initial.x, output.initial.y);

  // push initial gridcell onto queue
  std::queue<unsigned int> bfs;
  bfs.push(initial_cell);

  // cache reference position in world coords
  unsigned int rx, ry;
  double reference_x, reference_y;
  costmap_->indexToCells(reference, rx, ry);
  costmap_->mapToWorld(rx, ry, reference_x, reference_y);

  while (!bfs.empty()) {
    unsigned int idx = bfs.front();
    bfs.pop();

    // try adding cells in 8-connected neighborhood to frontier
    for (unsigned int nbr : nhood8(idx, *costmap_)) {
      // check if neighbour is a potential frontier cell
      if (isNewFrontierCell(nbr, frontier_flag)) {
        // mark cell as frontier
        frontier_flag[nbr] = true;
        unsigned int mx, my;
        double wx, wy;
        costmap_->indexToCells(nbr, mx, my);
        costmap_->mapToWorld(mx, my, wx, wy);

        geometry_msgs::Point point;
        point.x = wx;
        point.y = wy;
        output.points.push_back(point);

        // update frontier size
        output.size++;

        // update centroid of frontier
        output.centroid.x += wx;
        output.centroid.y += wy;

        // determine frontier's distance from robot, going by closest gridcell
        // to robot
        double distance = sqrt(pow((double(reference_x) - double(wx)), 2.0) +
                               pow((double(reference_y) - double(wy)), 2.0));
        if (distance < output.min_distance) {
          output.min_distance = distance;
          output.middle.x = wx;
          output.middle.y = wy;
        }

        // add to queue for breadth first search
        bfs.push(nbr);
      }
    }
  }

  // average out frontier centroid
  output.centroid.x /= output.size;
  output.centroid.y /= output.size;
  
  // [chad] draw eight points around the centroid and then pick any of them which are in free space
  //        the mean of those free-space points will be in the goal orientation
  // [chad] averaging angles may lead to ambiguous goals, so we average absolute positions and then calculate orientatioin
  double angles[8] = {0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0};
  double orientX[8], orientY[8];
  double orient_dist = 0.35;
  unsigned int orient_m_x, orient_m_y;
  int counted_vectors = 0;
  unsigned int orient_index;
  double goal_x = 0.0, goal_y = 0.0;
  double orient_x = 0.0, orient_y = 0.0;
  for (int i = 0; i < 8; i++){
      orientY[i] = output.centroid.y + orient_dist * sin(angles[i] * PI / 180.0);
      orientX[i] = output.centroid.x + orient_dist * cos(angles[i] * PI / 180.0);
      if (!costmap_->worldToMap(orientX[i], orientY[i], orient_m_x, orient_m_y)){
        ROS_DEBUG("[Warning] the orient is out of the map bounds %lf. Reducing orient_dist may prevent this warning", angles[i]);
	      counted_vectors += 1;
        goal_x += orientX[i];
        goal_y += orientY[i];
      }
      else{
	      orient_index = costmap_->getIndex(orient_m_x, orient_m_y);
        if (map_[orient_index] == NO_INFORMATION){
	        counted_vectors += 1;
	        goal_x += orientX[i];
	        goal_y += orientY[i];
        }
      }
  }
  if (counted_vectors != 0){
      goal_x = goal_x / counted_vectors;
      goal_y = goal_y / counted_vectors;
      orient_x = goal_x - output.centroid.x;
      orient_y = goal_y - output.centroid.y;
      if (orient_x != 0.0 || orient_y != 0.0){
	      output.theta = atan2(orient_y, orient_x);
      }
      else{
        ROS_DEBUG("[Warning] both arguments passed to atan2 are zero");
      }
  }
  else{
      ROS_DEBUG("[Warning] no open orients have been detected");
  }
  return output;
}

bool FrontierSearch::isNewFrontierCell(unsigned int idx,
                                       const std::vector<bool>& frontier_flag)
{
  // check that cell is unknown and not already marked as frontier
  if (map_[idx] != NO_INFORMATION || frontier_flag[idx]) {
    return false;
  }

  // frontier cells should have at least one cell in 4-connected neighbourhood
  // that is free
  for (unsigned int nbr : nhood4(idx, *costmap_)) {
    if (map_[nbr] == FREE_SPACE) {
      return true;
    }
  }

  return false;
}

double FrontierSearch::frontierCost(const Frontier& frontier)
{
  return (potential_scale_ * frontier.min_distance *
          costmap_->getResolution()) -
         (gain_scale_ * frontier.size * costmap_->getResolution());
}
}
