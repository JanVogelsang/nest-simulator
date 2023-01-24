#include "kernel_manager.h"


inline void
nest::SourceTablePosition::decrease()
{
  // first try finding a valid index by only decreasing target-local connection id
  --local_target_connection_id;
  if ( local_target_connection_id < 0 )
  {
    // then try finding a valid index by decreasing synapse index
    --syn_id;
    if ( syn_id < 0 )
    {
      // then try finding a valid index by decreasing target node id
      --local_target_node_id;
      if ( local_target_node_id < 0 )
      {
        // then try finding a valid index by decreasing thread index
        --tid;
        if ( tid < 0 )
        {
          return; // reached the end without finding a valid entry
        }

        local_target_node_id = kernel().node_manager.get_local_nodes( tid ).size() - 1;
      }

      syn_id = kernel().model_manager.get_num_connection_models() - 1;
    }

    local_target_connection_id = kernel()
                                   .node_manager.get_local_nodes( tid )
                                   .get_node_by_index( local_target_node_id )
                                   ->get_num_conn_type_sources( syn_id )
      - 1;
  }

  // if the found index is not valid, decrease further
  if ( local_target_connection_id
    == kernel()
         .node_manager.get_local_nodes( tid )
         .get_node_by_index( local_target_node_id )
         ->get_num_conn_type_sources( syn_id ) )
  {
    decrease();
  }
}

inline void
nest::SourceTablePosition::increase()
{
  // first try finding a valid index by only increasing target-local connection id
  ++local_target_connection_id;
  if ( local_target_connection_id
    == kernel()
         .node_manager.get_local_nodes( tid )
         .get_node_by_index( local_target_node_id )
         ->get_num_conn_type_sources( syn_id ) )
  {
    // then try finding a valid index by increasing synapse index
    ++syn_id;
    if ( syn_id == kernel().model_manager.get_num_connection_models() )
    {
      // then try finding a valid index by increasing target node id
      ++local_target_node_id;
      if ( local_target_node_id == kernel().node_manager.get_local_nodes( tid ).size() )
      {
        // then try finding a valid index by increasing thread index
        ++tid;
        if ( tid == kernel().vp_manager.get_num_threads() )
        {
          return; // reached the end without finding a valid entry
        }

        local_target_node_id = 0;
      }

      syn_id = 0;
    }

    local_target_connection_id = 0;
  }

  // if the found index is still not valid, increase further
  if ( local_target_connection_id
    == kernel()
         .node_manager.get_local_nodes( tid )
         .get_node_by_index( local_target_node_id )
         ->get_num_conn_type_sources( syn_id ) )
  {
    increase();
  }
}