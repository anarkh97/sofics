#include "RemoteSeedInitializerProxy.h"

namespace Pita {

RemoteSeedInitializerProxy::RemoteSeedInitializerProxy(Communicator * c, size_t vs) :
  SeedInitializer(vs),
  communicator_(c),
  sBuffer_(2 * vs),
  state_()
{}

DynamState
RemoteSeedInitializerProxy::initialSeed(SliceRank rank) const {
  StateMap::iterator it = state_.lower_bound(rank);
  if (it != state_.end() && it->first == rank) {
    return it->second;
  }

  communicator_->recFrom(rank.value(), sBuffer_.array(), 2 * vectorSize());
  DynamState newSeed = DynamState(vectorSize(), sBuffer_.array());
  state_.insert(it, std::make_pair(rank, newSeed));

  return newSeed;
}

} // end namespace Pita
