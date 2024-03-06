/*
 *  target_table_devices.cpp
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

// Includes from nestkernel:
#include "connection_label.h"
#include "connector_base.h"
#include "kernel_manager.h"
#include "node_impl.h"

namespace nest
{

TargetTableDevices::TargetTableDevices()
{
}

TargetTableDevices::~TargetTableDevices()
{
}

void
TargetTableDevices::initialize()
{
  const thread num_threads = kernel().vp_manager.get_num_threads();
  targets_to_devices_.resize( num_threads );
  targets_from_devices_.resize( num_threads );
  sending_devices_node_ids_.resize( num_threads );
}

void
TargetTableDevices::finalize()
{
  //  for ( size_t tid = 0; tid < targets_to_devices_.size(); ++tid )
  //  {
  //    for ( auto iit = targets_to_devices_[ tid ].begin(); iit != targets_to_devices_[ tid ].end(); ++iit )
  //    {
  //      for ( std::vector< ConnectorBase* >::iterator iiit = iit->begin(); iiit != iit->end(); ++iiit )
  //      {
  //        delete *iiit;
  //      }
  //    }
  //  }

  std::vector< std::map< index, std::vector< Target > > >().swap( targets_to_devices_ );
  std::vector< std::vector< std::vector< Target > > >().swap( targets_from_devices_ );
  std::vector< std::vector< index > >().swap( sending_devices_node_ids_ );
}

void
TargetTableDevices::resize_to_number_of_neurons()
{
#pragma omp parallel
  {
    const thread tid = kernel().vp_manager.get_thread_id();
    targets_from_devices_[ tid ].resize( kernel().node_manager.get_num_thread_local_devices( tid ) + 1 );
    sending_devices_node_ids_[ tid ].resize( kernel().node_manager.get_num_thread_local_devices( tid ) + 1 );
  } // end omp parallel
}

// void
// TargetTableDevices::resize_to_number_of_synapse_types()
//{
// #pragma omp parallel
//   {
//     const thread tid = kernel().vp_manager.get_thread_id();
//     //    for ( index lid = 0; lid < targets_to_devices_[ tid ].size(); ++lid )
//     //    {
//     //      // make sure this device has support for all synapse types
//     //      targets_to_devices_[ tid ][ lid ].resize( kernel().model_manager.get_num_connection_models(), nullptr );
//     //    }
//     for ( index ldid = 0; ldid < targets_from_devices_[ tid ].size(); ++ldid )
//     {
//       // make sure this device has support for all synapse types
//       targets_from_devices_[ tid ][ ldid ].resize( kernel().model_manager.get_num_connection_models() );
//     }
//   } // end omp parallel
// }

void
TargetTableDevices::get_synapse_status_from_device( const thread tid,
  const index ldid,
  const synindex syn_id,
  DictionaryDatum& dict,
  const index lcid ) const
{
  Node* target_node = kernel().node_manager.thread_lid_to_node( targets_from_devices_[ tid ][ ldid ][ lcid ].get_tid(),
    targets_from_devices_[ tid ][ ldid ][ lcid ].get_local_target_node_id() );
  target_node->get_connection_status( syn_id, lcid, dict );
}

void
TargetTableDevices::set_synapse_status_from_device( const thread tid,
  const index ldid,
  const synindex syn_id,
  ConnectorModel& cm,
  const DictionaryDatum& dict,
  const index lcid )
{
  Node* target_node = kernel().node_manager.thread_lid_to_node( targets_from_devices_[ tid ][ ldid ][ lcid ].get_tid(),
    targets_from_devices_[ tid ][ ldid ][ lcid ].get_local_target_node_id() );
  target_node->set_connection_status( syn_id, lcid, dict, cm );
}

}