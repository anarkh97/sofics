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


#include "allocators.h"
#include "contact_assert.h"
#include "ContactErrors.h"
#include "ContactFaceBlock.h"
#include "ContactNodeNodeInteraction.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactElementElementInteraction.h"
#include "ContactElement.h"
#include "ContactShellNode.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "contact_sorting.h"
#include "CString.h"
#include "ContactSearchData.h"
#include "search_methods.h"
#include "ContactParOStream.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactSequentialAllocator.h"
#include "ContactBoundingBox.h"
#include "ContactBoundingBoxHierarchy.h"
#include "ContactBoundingBoxHierarchy_Int.h"
#include "ContactRangeSearch.h"
#include "contact_tolerances.h"

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#include "ContactZoltan.h"
#endif

#include <iostream>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>

using namespace std;

void
ContactSearch::TrackedSearch(SearchType search_type, 
                             Real gap_tol, int num_configs,
                             VariableHandle CURRENT_POSITION, 
                             VariableHandle PREDICTED_POSITION,
                             VariableHandle AUGMENTED_POSITION )
{
  error_code = NO_ERROR;

  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( track_search_time );
  #endif
  
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
    if (global_tracking_interval>1) {
      postream<<"  Performing local tracked search, tracking_step = "<<tracking_step<<" of "<<global_tracking_interval-1<<"\n";
    } else {
      postream<<"  Performing local tracked search\n";
    }
    postream<<"    looking for "<<num_tracked_interactions
            <<" interactions on "<<num_tracked_nodes<<" nodes\n";
  #endif
  #ifdef CONTACT_HEARTBEAT 
    int num_tracked_nodes_global        = contact_global_sum(num_tracked_nodes,  SearchComm);
    int num_tracked_interactions_global = contact_global_sum(num_tracked_interactions, SearchComm);
    if (contact_processor_number(SearchComm)==0) {
      if (global_tracking_interval>1) {
        std::cout<<"  Performing local tracked search, tracking_step = "<<tracking_step<<" of "<<global_tracking_interval-1<<"\n";
      } else {
        std::cout<<"  Performing local tracked search\n";
      }
      std::cout<<"    looking for "<<num_tracked_interactions_global
               <<" interactions on "<<num_tracked_nodes_global<<" nodes\n";
    }
  #endif
        
  // Tracked (local) search is done in the primary decomposition
  search_topology = primary_topology;

  #ifndef CONTACT_NO_MPI
  if( contact_number_of_processors(SearchComm)>1) {
    #ifdef CONTACT_TIMINGS
      timer.Start_Timer( track_ghost_time );
    #endif
    if (keep_ghosting) {
      if (tracking_step==1) {
        #ifdef CONTACT_TIMINGS
          timer.Start_Timer( track_cleanup_time );
        #endif
        #if CONTACT_DEBUG_PRINT_LEVEL>=4
          postream<<"    calling DeleteGhostedNFIfaces()\n";
        #else
        #ifdef CONTACT_HEARTBEAT 
          if (contact_processor_number(SearchComm)==0) {
            std::cout<<"    calling DeleteGhostedNFIfaces()\n";
          }
        #endif
        #endif
        primary_topology->DeleteGhostedNFIfaces();
        #ifdef CONTACT_TIMINGS
          timer.Stop_Timer( track_cleanup_time );
        #endif
        int flag = global_tracking_interval>1?1:0;
        #if CONTACT_DEBUG_PRINT_LEVEL>=4
          postream<<"    calling GhostNFIfaces("<<flag<<")\n";
        #else
        #ifdef CONTACT_HEARTBEAT 
          if (contact_processor_number(SearchComm)==0) {
            std::cout<<"    calling GhostNFIfaces("<<flag<<")\n";
          }
        #endif
        #endif
	error_code = primary_topology->GhostNFIfaces(flag);
	if (contact_global_error_check(error_code, SearchComm)) {
	  errors->Add_Error_Message("Error in ghosting, possibly due to");
	  errors->Add_Error_Message("a failure to re-connect a constraint");
	  errors->Add_Error_Message("node and face where one is ghosted");
	  return;
	}
      } else {
        #if CONTACT_DEBUG_PRINT_LEVEL>=4
          postream<<"    calling UpdateGhostedNFIfaces()\n";
        #else
        #ifdef CONTACT_HEARTBEAT 
          if (contact_processor_number(SearchComm)==0) {
            std::cout<<"    calling UpdateGhostedNFIfaces()\n";
          }
        #endif
        #endif
        primary_topology->UpdateGhostedNFIfaces();
      }
    } else {
      int flag = 0;
      #if CONTACT_DEBUG_PRINT_LEVEL>=4
        postream<<"    calling GhostNFIfaces("<<flag<<")\n";
      #endif
      #ifdef CONTACT_HEARTBEAT 
        if (contact_processor_number(SearchComm)==0) {
          std::cout<<"    calling GhostNFIfaces("<<flag<<")\n";
        }
      #endif
      error_code = primary_topology->GhostNFIfaces(flag);
      if (contact_global_error_check(error_code, SearchComm)) {
	errors->Add_Error_Message("Error in ghosting, possibly due to");
	errors->Add_Error_Message("a failure to re-connect a constraint");
	errors->Add_Error_Message("node and face where one is ghosted");
	return;
      }
    }
    #ifdef CONTACT_TIMINGS
      timer.Stop_Timer( track_ghost_time );
    #endif

  }
  #endif
  
  // Retrieve any tied interactions
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( search_retr_tied_time );
  #endif
  switch (num_configs) {
  case 1:
    Retrieve_Tied_Interactions_From_Primary(CURRENT_POSITION);
    break;
  case 2:
    Retrieve_Tied_Interactions_From_Primary(PREDICTED_POSITION);
    break;
  case 3:
    Retrieve_Tied_Interactions_From_Primary(PREDICTED_POSITION);
    break;
  default:
    POSTCONDITION(0);
    break;
  }
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( search_retr_tied_time );
  #endif

  // Retrieve any glued interactions
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( search_retr_tied_time );
  #endif
  switch (num_configs) {
  case 1:
      Retrieve_Glued_Interactions(CURRENT_POSITION, ContactSearch::SEARCH, 1);
      break;
  case 2:
      Retrieve_Glued_Interactions(PREDICTED_POSITION, ContactSearch::SEARCH, 1);
      break;
  case 3:
      Retrieve_Glued_Interactions(PREDICTED_POSITION, ContactSearch::SEARCH, 1);
      break;
  default:
      POSTCONDITION(0);
  break;
  }
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( search_retr_tied_time );
  #endif

  // Retrieve any tracked interactions
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( track_retr_tracked_time );
  #endif
  Retrieve_Tracked_Interactions(search_type, gap_tol);
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( track_retr_tracked_time );
  #endif

  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( track_id_time );
  #endif
  Interaction_Definition( num_configs, ContactSearch::PRIMARY, 
                          ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE );
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( track_id_time );
  #endif
      
  #ifndef CONTACT_NO_MPI 
  if( contact_number_of_processors(SearchComm)>1) {
    if (!keep_ghosting) {
      #ifdef CONTACT_TIMINGS
        timer.Start_Timer( track_cleanup_time );
      #endif
      #if CONTACT_DEBUG_PRINT_LEVEL>=4
        postream<<"    DeleteGhostedNFIfaces()\n";
      #else
      #ifdef CONTACT_HEARTBEAT 
        if (contact_processor_number(SearchComm)==0) {
          std::cout<<"    DeleteGhostedNFIfaces()\n";
        }
      #endif
      #endif
      primary_topology->DeleteGhostedNFIfaces();
      #ifdef CONTACT_TIMINGS
        timer.Stop_Timer( track_cleanup_time );
      #endif
    }
  }
  #endif

#if CONTACT_DEBUG_PRINT_LEVEL>=2
  primary_topology->Display_NodeEntity_Interactions_Summary(
      postream, ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE, (char*)"    ", 0);
#else
#ifdef CONTACT_HEARTBEAT       
  primary_topology->Display0_NodeEntity_Interactions_Summary(
      ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE, (char*)"    ", 0 );
#endif
#endif

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( track_search_time );
#endif
}
