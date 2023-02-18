#include "adjacency_list.h"

namespace nest
{

void AdjacencyList::add_target( const thread tid, const index source_node_id, const index target_node_id, const index target_connection_id )
{
  auto index = sources_[ tid ].find( source_node_id );

  if ( index != sources_[ tid ].end() )
  {
    adjacency_list_[ tid ][ (*index).second ].emplace_back( target_node_id, target_connection_id );
  }
  else {
    sources_[ tid ][ source_node_id ] = adjacency_list_.size();
    adjacency_list_[ tid ].emplace_back( std::initializer_list< AdjacencyListTarget >{ { target_node_id, target_connection_id } } );
  }
}

} // nest