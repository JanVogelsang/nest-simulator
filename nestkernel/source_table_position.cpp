/*
 *  source_table_position.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "source_table_position.h"

// Includes from nestkernel:
#include "kernel_manager.h"
#include "mpi_manager_impl.h"


namespace nest
{

void
SourceTablePosition::decrease( const thread source_rank )
{
  // we stay in this loop either until we can found a valid target or we have reached the end of all sources tables
  do
  {
    // first try finding a valid index by only decreasing target-local connection id
    --local_target_connection_id;

    // keep decreasing until a valid index is found
    while ( local_target_connection_id < 0 )
    {
      // then try finding a valid index by decreasing synapse index
      --syn_id;
      if ( syn_id < 0 )
      {
        Node* target_node;
        // then try finding a valid index by decreasing target node id
        do
        {
          --local_target_node_id;

          while ( local_target_node_id < 0 )
          {
            // then try finding a valid index by decreasing thread index
            --tid;
            if ( tid < 0 )
            {
              tid = -1;
              local_target_node_id = -1;
              syn_id = -1;
              local_target_connection_id = -1;
              return; // reached the end without finding a valid entry
            }

            local_target_node_id = kernel().node_manager.get_local_nodes( tid ).size() - 1;
          }

          target_node = kernel().node_manager.get_local_nodes( tid ).get_node_by_index( local_target_node_id );
        } while ( not target_node->has_proxies() or target_node->get_node_id() == DISABLED_NODE_ID ); // skip devices

        syn_id = kernel().model_manager.get_num_connection_models() - 1;
      }

      local_target_connection_id = kernel()
                                     .node_manager.get_local_nodes( tid )
                                     .get_node_by_index( local_target_node_id )
                                     ->get_num_conn_type_sources( syn_id )
        - 1;
    }
    current_source = kernel()
                       .node_manager.get_local_nodes( tid )
                       .get_node_by_index( local_target_node_id )
                       ->get_source( syn_id, local_target_connection_id );

    thread current_source_rank = kernel().mpi_manager.get_process_id_of_node_id( current_source );
    if ( source_rank == current_source_rank and current_source != DISABLED_NODE_ID )
    {
      return; // found a valid entry
    }
  } while ( not reached_end() );
}

void
SourceTablePosition::increase( const thread source_rank )
{
  do
  {
    // first try finding a valid index by only increasing target-local connection id
    ++local_target_connection_id;

    // keep increasing until a valid index is found
    while ( local_target_connection_id
      == static_cast< long >( kernel()
                                .node_manager.get_local_nodes( tid )
                                .get_node_by_index( local_target_node_id )
                                ->get_num_conn_type_sources( syn_id ) ) )
    {
      if ( local_target_connection_id
        == static_cast< long >( kernel()
                                  .node_manager.get_local_nodes( tid )
                                  .get_node_by_index( local_target_node_id )
                                  ->get_num_conn_type_sources( syn_id ) ) )
      {
        // then try finding a valid index by increasing synapse index
        ++syn_id;
        if ( syn_id == static_cast< long >( kernel().model_manager.get_num_connection_models() ) )
        {
          Node* target_node;
          // then try finding a valid index by increasing target node id
          do
          {
            ++local_target_node_id;

            while ( local_target_node_id == static_cast< long >( kernel().node_manager.get_local_nodes( tid ).size() ) )
            {
              // then try finding a valid index by increasing thread index
              ++tid;
              if ( tid == static_cast< long >( kernel().vp_manager.get_num_threads() ) )
              {
                tid = -1;
                local_target_node_id = -1;
                syn_id = -1;
                local_target_connection_id = -1;
                return; // reached the end without finding a valid entry
              }

              local_target_node_id = 0;
            }

            target_node = kernel().node_manager.get_local_nodes( tid ).get_node_by_index( local_target_node_id );
          } while ( not target_node->has_proxies() or target_node->get_node_id() == DISABLED_NODE_ID ); // skip devices

          syn_id = 0;
        }

        local_target_connection_id = 0;
      }
    }
    // the current position contains an entry, so we retrieve it
    current_source = kernel()
                       .node_manager.get_local_nodes( tid )
                       .get_node_by_index( local_target_node_id )
                       ->get_source( syn_id, local_target_connection_id );

    thread current_source_rank = kernel().mpi_manager.get_process_id_of_node_id( current_source );
    if ( source_rank == current_source_rank and current_source != DISABLED_NODE_ID )
    {
      return; // found a valid entry
    }
  } while ( not reached_end() );
}

}