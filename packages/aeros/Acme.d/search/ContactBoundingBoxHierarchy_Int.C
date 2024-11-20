// This file is part of a modified version of ACME: Algorithms for Contact in
// a Multiphysics Environment, derived from version 2.7f
//
// Copyright (C) 2007 Sandia Corporation
// Copyright (C) 2011 Stanford University 
//
// ACME is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ACME.  If not, see <http://www.gnu.org/licenses/>.


#include <ContactBoundingBoxHierarchy_Int.h>
#include <Contact_Communication.h>
#include "ContactTimer.h"
#include "CString.h"
#include "contact_sorting.h"
#include "cmath"

using namespace std;

int*                          ObjectBoundingBoxHierarchy_Int::index_array = NULL;
int*                          ObjectBoundingBoxHierarchy_Int::index_array0 = NULL;
int*                          ObjectBoundingBoxHierarchy_Int::index_array1 = NULL;
int*                          ObjectBoundingBoxHierarchy_Int::index_array2 = NULL;

static inline void split_boxes(const int &num_boxes, const ObjectBoundingBox_Int *boxes, float *moment) {
  static float global_center[3];
  static float centroid[3];
  //
  //  Calculate the centroid of each component box, the centroid of the current encompassing box, and the 
  //  moment of inertia of all particles in the box about (0, 0, 0)
  //
  global_center[0] = 0.0;
  global_center[1] = 0.0;
  global_center[2] = 0.0;
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    boxes[ibox].calculate_centroid2(centroid);
    global_center[0] += centroid[0];
    global_center[1] += centroid[1];
    global_center[2] += centroid[2];
    moment[0] += centroid[0] * centroid[0];
    moment[1] += centroid[1] * centroid[1];
    moment[2] += centroid[2] * centroid[2];
  } 
  //
  //  Adjust the moment to be taken about the global_center location
  //  Note:  Here we are actually calculating num_boxes * moment.  Only need to compare relative values
  //  of moment so multiplying by num_boxes is fine.  This is done to avoid any division operations.  
  //
  moment[0] *= num_boxes;
  moment[0] -= global_center[0] * global_center[0];
  moment[1] *= num_boxes;
  moment[1] -= global_center[1] * global_center[1];
  moment[2] *= num_boxes;
  moment[2] -= global_center[2] * global_center[2];
}


ObjectBoundingBoxHierarchy_Int::ObjectBoundingBoxHierarchy_Int() :
right_child_offset(0)
{
}

ObjectBoundingBoxHierarchy_Int::~ObjectBoundingBoxHierarchy_Int() {
}

void ObjectBoundingBoxHierarchy_Int::create_hierarchy(ObjectBoundingBoxHierarchy_Int *hierarchy_data, 
                                                  ObjectBoundingBox_Int *boxes, 
                                                  const int &num_boxes) {
  //
  //  Allocate some data arrays nessecary for creating the bounding box hierarchy.  Arrays are alloated to allow creating a MAX_TREE_LEVELS
  //  level hierarchy, this should be sufficently large to handle in concievable problem.
  //
  index_array          = new int[num_boxes];
  index_array0         = new int[num_boxes];
  index_array1         = new int[num_boxes];
  index_array2         = new int[num_boxes];

  float *cent_array = new float[num_boxes * 3];

  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    float centroid[3];
    boxes[ibox].calculate_centroid2(centroid);
    cent_array[ibox]               = centroid[0];
    cent_array[ibox + num_boxes]   = centroid[1];
    cent_array[ibox + num_boxes*2] = centroid[2];
  }
  FORTRAN(contact_indexx_float)(num_boxes, cent_array,                 index_array0, num_boxes);
  FORTRAN(contact_indexx_float)(num_boxes, &(cent_array[num_boxes]),   index_array1, num_boxes);
  FORTRAN(contact_indexx_float)(num_boxes, &(cent_array[num_boxes*2]), index_array2, num_boxes);
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    index_array0[ibox]--;
    index_array1[ibox]--;
    index_array2[ibox]--;
  }

  //
  //  Create a list of the boxes that hierarchy will contain
  //

  ObjectBoundingBox_Int *scratch_boxes = new ObjectBoundingBox_Int[num_boxes];
  hierarchy_data->create_hierarchy_loop(hierarchy_data,
                                        boxes, 
                                        scratch_boxes, 
                                        num_boxes,
                                        index_array0,
                                        index_array1,
                                        index_array2);
  delete [] cent_array;
  delete [] scratch_boxes;
  //
  //  Delete all static data items
  //
  delete [] index_array;
  delete [] index_array0;
  delete [] index_array1;
  delete [] index_array2;
}

struct CreateHierarchy_IntStack {
  ObjectBoundingBoxHierarchy_Int *hierarchy_loc;
  ObjectBoundingBox_Int *boxes;
  ObjectBoundingBox_Int *scratch_boxes;
  int num_boxes;
  int *index_array0_t;
  int *index_array1_t;
  int *index_array2_t;
  int unwind;
};

//
//  Store a one noded tree, just set the terminal case
//
static inline void store_1_node_tree(ObjectBoundingBoxHierarchy_Int *root_node, const ObjectBoundingBox_Int *const boxes) {
  root_node->set_right_child_offset(-boxes[0].get_object_number());
  root_node->set_box(boxes[0]);
}

static inline void store_2_node_tree(ObjectBoundingBoxHierarchy_Int *root_node, const ObjectBoundingBox_Int *const boxes) {
  //
  //  Special case, only two boxes passed in.  This happens often due to the top to bottom tree structure.
  //  Significant optimizations are avalable for this case
  //  Split the list into two halves.  Put half the boxes into one child and half into the other.  Complete processing
  //  on both halves.
  //
  //  Store self
  //
  root_node->set_right_child_offset(2);
  root_node->set_box(boxes[0]);
  root_node->add_box(boxes[1]);
  //
  //  Store left child
  //
  root_node[1].set_right_child_offset(-boxes[0].get_object_number());
  root_node[1].set_box(boxes[0]);
  //
  //  Store right child
  //
  root_node[2].set_right_child_offset(-boxes[1].get_object_number());
  root_node[2].set_box(boxes[1]);
}

static inline void store_3_node_tree(ObjectBoundingBoxHierarchy_Int *root_node, const ObjectBoundingBox_Int *const boxes) {
  //
  //  Special case, exatly three boxes pased in.  This happens resonably often, aditionally, if handle 2 box case, and
  //  3 box case then automatically will handle the one box cases (except for the very special case of exactly one box
  //  in the entire hierarchy, which is handled in the calling routine.)
  //
  //  Determine how to split the boxes.  The closest two boxes will be placed in the right child,
  //  and the third box in the left child.  The index arrays will not be used for this operation
  //  for optimization purposes
  //
  float dist_squared_01 = abs((boxes[0].get_x_min() + boxes[0].get_x_max()) - (boxes[1].get_x_min() + boxes[1].get_x_max()) + 
			      (boxes[0].get_y_min() + boxes[0].get_y_max()) - (boxes[1].get_y_min() + boxes[1].get_y_max()) +
			      (boxes[0].get_z_min() + boxes[0].get_z_min()) - (boxes[1].get_z_min() + boxes[1].get_z_min()));
  float dist_squared_12 = abs((boxes[1].get_x_min() + boxes[1].get_x_max()) - (boxes[2].get_x_min() + boxes[2].get_x_max()) + 
			      (boxes[1].get_y_min() + boxes[1].get_y_max()) - (boxes[2].get_y_min() + boxes[2].get_y_max()) +
			      (boxes[1].get_z_min() + boxes[1].get_z_min()) - (boxes[2].get_z_min() + boxes[2].get_z_min()));
  float dist_squared_20 = abs((boxes[2].get_x_min() + boxes[2].get_x_max()) - (boxes[0].get_x_min() + boxes[0].get_x_max()) + 
			      (boxes[2].get_y_min() + boxes[2].get_y_max()) - (boxes[0].get_y_min() + boxes[0].get_y_max()) +
			      (boxes[2].get_z_min() + boxes[2].get_z_min()) - (boxes[0].get_z_min() + boxes[0].get_z_min()));
  //
  //  Create the left and right hierarchies
  //
  if(dist_squared_01 < dist_squared_12 && dist_squared_01 < dist_squared_20) {
    //
    //  Right child has boxes 0 and 1, left child has box 2
    //
    (root_node+1)->set_right_child_offset(-boxes[2].get_object_number());
    (root_node+1)->set_box(boxes[2]);
    (root_node+3)->set_right_child_offset(-boxes[0].get_object_number());
    (root_node+3)->set_box(boxes[0]);
    (root_node+4)->set_right_child_offset(-boxes[1].get_object_number());
    (root_node+4)->set_box(boxes[1]);
    (root_node+2)->set_right_child_offset(2);      
    (root_node+2)->set_box(boxes[0]);
    (root_node+2)->add_box(boxes[1]);  
  } else {
    if(dist_squared_12 < dist_squared_20) {
      //
      //  Right child has boxes 1 and 2, left child has box 0
      //
      (root_node+1)->set_right_child_offset(-boxes[0].get_object_number());
      (root_node+1)->set_box(boxes[0]);
      (root_node+3)->set_right_child_offset(-boxes[1].get_object_number());
      (root_node+3)->set_box(boxes[1]);
      (root_node+4)->set_right_child_offset(-boxes[2].get_object_number());
      (root_node+4)->set_box(boxes[2]);
      (root_node+2)->set_right_child_offset(2);      
      (root_node+2)->set_box(boxes[1]);
      (root_node+2)->add_box(boxes[2]);  
    } else {
      //
      //  Right child has boxes 0 and 2, left child has box 1
      //
      (root_node+1)->set_right_child_offset(-boxes[1].get_object_number());
      (root_node+1)->set_box(boxes[1]);
      (root_node+3)->set_right_child_offset(-boxes[0].get_object_number());
      (root_node+3)->set_box(boxes[0]);
      (root_node+4)->set_right_child_offset(-boxes[2].get_object_number());
      (root_node+4)->set_box(boxes[2]);
      (root_node+2)->set_right_child_offset(2);
      (root_node+2)->set_box(boxes[0]);
      (root_node+2)->add_box(boxes[2]);  
    }
  }
  root_node->set_right_child_offset(2);      
  root_node->set_box(*(root_node+1));
  root_node->add_box(*(root_node+2));
}

//
//  Non recursive hierarchy creation routine
//
//  Removing recurse components improves speed.  However, recursion is easier to understand, here is the 
//  Recursive puesdo code for this routine:
//
//  ObjectBoundingBoxHierarchy_Int::create_hierarchy(boxes, scratch_boxes, num_boxes, index0, index1, index2) {
//    if(num_boxes == 1) {
//      1) Set the bounding box of this tree node to the bounding box of the object that the tree node represents.
//      2) Set the right child offset of this node to the negative of the object number this node represents
//    }
//    if(num_boxes > 1) {
//      1) Calculate the number of boxes to place into the left and right children
//      2) Determine which direction to split these boxes, this direction is taken as the direction with
//         the largest moment of inertia of box centroids
//      3) Sort the boxes based on the split direction using the index array for that direction.  Place sorted
//         boxes into scratch_boxes.
//      4) Update the sorting lists and create the left and right child sorting lists.
//      5) Recurively call the routine for the left and right child
//      this[1].create_hierachy(left_boxes, left_scratch_boxes, left_num_boxes, left_index0, left_index1, left_index2)
//      this[right_child_offset].create_hierachy(right_boxes, right_scratch_boxes, right_num_boxes, right_index0, right_index1, right_index2)
//      6) Recursive unwinding, set the bounding box of the node as the sum of the bounding boxes of the child nodes
//    }
//  }
//
void ObjectBoundingBoxHierarchy_Int::create_hierarchy_loop(ObjectBoundingBoxHierarchy_Int *const hierarchy_start_ptr,
                                                       ObjectBoundingBox_Int *const boxes_start_ptr, 
                                                       ObjectBoundingBox_Int *const scratch_boxes_start_ptr,
                                                       const int num_boxes_start,
                                                       int *const index_array0_t_start_ptr,
                                                       int *const index_array1_t_start_ptr,
                                                       int *const index_array2_t_start_ptr) {
  //
  //  Create the object stack variables, this routine is really a recursive routine, however, some performance gain
  //  can be had by rewriting as a loop and managing the recursive stack by hand.
  //
  CreateHierarchy_IntStack object_stack[MAX_TREE_LEVELS];
  CreateHierarchy_IntStack *stack_ptr = object_stack;
  CreateHierarchy_IntStack *current_object = stack_ptr;
  //
  //  The tree creation loop handles several special cases, the effect is that the case of one, two and three boxs is not handled.  
  //  correctly, if it occurs handle it here and exit tree creation
  //
  if(num_boxes_start == 1) {
    store_1_node_tree(hierarchy_start_ptr, boxes_start_ptr);
    return;
  } else if(num_boxes_start ==2) {
    store_2_node_tree(hierarchy_start_ptr, boxes_start_ptr);
    return;
  } else if(num_boxes_start ==3) {
    store_3_node_tree(hierarchy_start_ptr, boxes_start_ptr);
    return;
  }
  //
  //  Push the head object onto the stack, this object contains all input boxes and will represent the root
  //  node of the object tree.
  //
  current_object->hierarchy_loc = hierarchy_start_ptr;
  current_object->boxes = boxes_start_ptr;
  current_object->scratch_boxes = scratch_boxes_start_ptr;
  current_object->num_boxes = num_boxes_start;
  current_object->index_array0_t = index_array0_t_start_ptr;
  current_object->index_array1_t = index_array1_t_start_ptr;
  current_object->index_array2_t = index_array2_t_start_ptr;
  current_object->unwind = 0;

  do {
    if(current_object->unwind) {
      //
      //  Stack unwinding step.  Compute the current node bounding box from the sum of the two child boxes.  Note, the last step
      //  of tree creation will be the unwinding of the root node, thus check for the terminal condition here.
      //
      ObjectBoundingBoxHierarchy_Int *hierarchy_loc = current_object->hierarchy_loc;      
      hierarchy_loc->set_box(*(hierarchy_loc+1));
      hierarchy_loc->add_box(*(hierarchy_loc+hierarchy_loc->right_child_offset));
      if(stack_ptr != object_stack) {
        current_object = --stack_ptr;
        continue;
      } else {
        break;
      }
    } else {
      const int num_boxes = current_object->num_boxes;     
      if(num_boxes ==2) {
        store_2_node_tree(current_object->hierarchy_loc, current_object->boxes);
        current_object = --stack_ptr;
        continue;
      } else if(num_boxes == 3) {
        store_3_node_tree(current_object->hierarchy_loc, current_object->boxes);
        current_object = --stack_ptr;
        continue;
      } else if (num_boxes <= 6){ 
	//
	//  Special case, do not update splitting directions, the unique cases for two and three boxes do not use them
	//  and splitting a N<=6 boxes yields only 3, or 2 length sub trees.
	//
	ObjectBoundingBoxHierarchy_Int *hierarchy_loc = current_object->hierarchy_loc;      
	int *const index_array0_t = current_object->index_array0_t;
	int *const index_array1_t = current_object->index_array1_t;
	int *const index_array2_t = current_object->index_array2_t;
	ObjectBoundingBox_Int *const boxes = current_object->boxes; 
	ObjectBoundingBox_Int *const scratch_boxes = current_object->scratch_boxes;
	//
	//  There are more than 2 boxes, compute an optimal splitting direction for the boxes, and divide them into two halves.
	//  Compute the centroid of each component box, the centroid of the current encompasing box, and the moment of 
	//  interial of all particles in the box about (0,0,0)
	//
	const int right_child_size = num_boxes/2;
	const int left_child_size = num_boxes - right_child_size;
	float moment[3]        = {0.0, 0.0, 0.0};
	split_boxes(num_boxes, boxes, moment);
	//
	//  Determine the longest centroid bounding box direction.  This is the direction in which bounding boxes will be 
	//  sorted and split.  Reorder the box sorting arrays based on the split direction
	//
	int *master_index;
	if(moment[0] > moment[1] && moment[0] > moment[2]) {
	  master_index = index_array0_t;
	} else {
	  if(moment[1] > moment[2]) {
	    master_index = index_array1_t;
	  } else {
	    master_index = index_array2_t;
	  }
	}
	//
	//  Sort boxes into array scratch_boxes
	//  Create the new master index array
	//
	for(int ibox = 0; ibox < num_boxes; ++ibox) {
	  scratch_boxes[ibox] = boxes[master_index[ibox]];
	}
	//
	//  Update the current tree pointer for right child offset
	//
	const int right_child_offset = left_child_size*2;
	hierarchy_loc->right_child_offset = right_child_offset;
	//
	//  Primary computations on this object are complete, thus set the unwinding flag
	//
	current_object->unwind = 1;
	//
	//  Add the right child to the object stack
	//
	current_object = ++stack_ptr;      
	current_object->hierarchy_loc = hierarchy_loc+right_child_offset;
	current_object->boxes = scratch_boxes + left_child_size;
	current_object->num_boxes = right_child_size;
	current_object->unwind = 0;
	//
	//  Add the left child box to the stack,
	// 
	current_object = ++stack_ptr;
	//
	//  Update the current object with the rest of the left objects data
	//
	current_object->hierarchy_loc = hierarchy_loc+1;
	current_object->unwind = 0;
	current_object->boxes = scratch_boxes;
	current_object->num_boxes = left_child_size;
	current_object->unwind = 0;
	continue;
      } else {
	ObjectBoundingBoxHierarchy_Int *hierarchy_loc = current_object->hierarchy_loc;      
	int *const index_array0_t = current_object->index_array0_t;
	int *const index_array1_t = current_object->index_array1_t;
	int *const index_array2_t = current_object->index_array2_t;
	ObjectBoundingBox_Int *const boxes = current_object->boxes; 
	ObjectBoundingBox_Int *const scratch_boxes = current_object->scratch_boxes;

	//
	//  There are more than 6 boxes, compute an optimal splitting direction for the boxes, and divide them into two halves.
	//  Compute the centroid of each component box, the centroid of the current encompasing box, and the moment of 
	//  interial of all particles in the box about (0,0,0)
	//
	const int right_child_size = num_boxes/2;
	const int left_child_size = num_boxes - right_child_size;
	float moment[3]        = {0.0, 0.0, 0.0};
	split_boxes(num_boxes, boxes, moment);
	//
	//  Determine the longest centroid bounding box direction.  This is the direction in which bounding boxes will be 
	//  sorted and split.  Reorder the box sorting arrays based on the split direction
	//
	int *master_index, *sub_index1, *sub_index2;
	if(moment[0] > moment[1] && moment[0] > moment[2]) {
	  master_index = index_array0_t;
	  sub_index1   = index_array1_t;
	  sub_index2   = index_array2_t;
	} else {
	  if(moment[1] > moment[2]) {
	    master_index = index_array1_t;
	    sub_index1   = index_array0_t;
	    sub_index2   = index_array2_t;
	  } else {
	    master_index = index_array2_t; 
	    sub_index1   = index_array0_t;
	    sub_index2   = index_array1_t;
	  }
	}
	//
	//  Sort boxes into array scratch_boxes
	//  Create the new master index array
	//
	for(int ibox = 0; ibox < num_boxes; ++ibox) {
	  scratch_boxes[ibox] = boxes[master_index[ibox]];
	  index_array[master_index[ibox]] = ibox; 
	}
	//
	//  Reorder secondary arrays to be consistent with the master array splitting
	//
	int left_pos1 = 0;
	int right_pos1 = 0; 
	int left_pos2 = 0;
	int right_pos2 = 0; 

	int *master_temp = master_index + right_child_size;

	for(int ibox = 0; ibox < num_boxes; ++ibox) {
	  const int master1 = index_array[sub_index1[ibox]];
	  const int master2 = index_array[sub_index2[ibox]];

	  if(master1 < left_child_size) {
	    sub_index1[left_pos1++] = master1;
	  } else {
	    master_index[right_pos1++] = master1 - left_child_size;
	  }

	  if(master2 < left_child_size) {
	    sub_index2[left_pos2++] = master2;
	  } else {
	    master_temp[right_pos2++] = master2 - left_child_size;
	  }
	}
	for(int ibox = 0; ibox < right_child_size; ++ibox) {
	  sub_index1[left_child_size + ibox] = master_index[ibox];
	  sub_index2[left_child_size + ibox] = master_temp[ibox];
	}
	//
	//  Update the master index to be correct with the new box ordering
	//
	for(int ibox = 0; ibox < left_child_size; ++ibox) {
	  master_index[ibox] = ibox;
	}
	for(int ibox = left_child_size; ibox < num_boxes; ++ibox) {
	  master_index[ibox] = ibox - left_child_size;
	}
	//
	//  Update the current tree pointer for right child offset
	//
	const int right_child_offset = left_child_size*2;
	hierarchy_loc->right_child_offset = right_child_offset;
	//
	//  Primary computations on this object are complete, thus set the unwinding flag
	//
	current_object->unwind = 1;
	//
	//  Add the right child to the object stack
	//
	current_object = ++stack_ptr;      
	current_object->hierarchy_loc = hierarchy_loc+right_child_offset;
	current_object->boxes = scratch_boxes + left_child_size;
	current_object->scratch_boxes = boxes + left_child_size;
	current_object->num_boxes = right_child_size;
	current_object->index_array0_t = index_array0_t + left_child_size;
	current_object->index_array1_t = index_array1_t + left_child_size;
	current_object->index_array2_t = index_array2_t + left_child_size;
	current_object->unwind = 0;
	//
	//  Add the left child box to the stack,
	// 
	current_object = ++stack_ptr;
	//
	//  Update the current object with the rest of the left objects data
	//
	current_object->hierarchy_loc = hierarchy_loc+1;
	current_object->unwind = 0;
	current_object->boxes = scratch_boxes;
	current_object->scratch_boxes = boxes;
	current_object->num_boxes = left_child_size;
	current_object->index_array0_t = index_array0_t;
	current_object->index_array1_t = index_array1_t;
	current_object->index_array2_t = index_array2_t;
	current_object->unwind = 0;
	continue;
      }
    }
  } while(true); 
}


//
//  NKC, This is a non recursive version of "search_for_overlap_recurse", Recursion removed in an attempt to 
//  Elminate the function overhead and operation system stack overhead.
//
//  Recursive version easier to understand, recurisve psuedo code:
//
//  ObjectBoundingBoxHierarchy_Int::search(object_box, obj_list, num_objs) {
//    if object box does not overlap the current object, return
//    if terminal node (right child offset < 1), add current node to the object list and return
//    if not terminal node, recurse to child nodes
//      this[1].search(object_box, obj_list, num_objs);
//      this[right_child_offset].search(object_box, obj_list, num_objs);
//  }
//

void ObjectBoundingBoxHierarchy_Int::search_for_overlap_loop(const ObjectBoundingBoxHierarchy_Int *const hierarchy_start_ptr,
							 const ContactBoundingBox_Int &object_box, 
							 int *obj_list, 
							 int &num_objs) {
  ObjectBoundingBoxHierarchy_Int const* current_object  = hierarchy_start_ptr;
  ObjectBoundingBoxHierarchy_Int const* object_stack[MAX_TREE_LEVELS];
  const ObjectBoundingBoxHierarchy_Int ** stack_ptr = object_stack;

  do {
    if(current_object->overlap(object_box)) {
      const int right_child_offset = current_object->right_child_offset;
      if(right_child_offset > 0) {
        *(stack_ptr++) = current_object++ + right_child_offset;
        continue;
      }
      obj_list[num_objs++] = -right_child_offset;
    }
    if(stack_ptr != object_stack) {
      current_object = *(--stack_ptr);
      continue;
    } 
    break;
  } while(true);
}



//
//  Search the box hierarchy.  Determine which objects the object overlaps.  Place the overlaping objects in a list.
//  If node is not a terminal node recursivley search all nodes owned by the current node.
//
void ObjectBoundingBoxHierarchy_Int::search_for_overlap_recurse_sym(const ObjectBoundingBoxHierarchy_Int &object_box, 
								int *obj_list, 
								int &num_objs,
								const ObjectBoundingBoxHierarchy_Int *max_owned_ptr) {
  if(max_owned_ptr <= &object_box) return;
  if(!overlap(object_box)) return;
  if(right_child_offset <= 0) { 
    obj_list[num_objs] = -right_child_offset;
    num_objs++;
  } else {
    this[1].search_for_overlap_recurse_sym(object_box, obj_list, num_objs, &(this[right_child_offset -1]));
    this[right_child_offset].search_for_overlap_recurse_sym(object_box, obj_list, num_objs, max_owned_ptr);
  }
}

