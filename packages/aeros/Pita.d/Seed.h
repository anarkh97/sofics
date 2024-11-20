#ifndef PITA_SEED_H
#define PITA_SEED_H

#include "SharedState.h"

#include "DynamState.h"
#include <Math.d/Vector.h>

namespace Pita {

typedef SharedState<DynamState> Seed;
typedef SharedState<Vector> ReducedSeed;

} // end namespace Pita

#endif /* PITA_SEED_H */
