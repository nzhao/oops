#ifndef OOPS_H
#define OOPS_H

#ifdef HAS_MATLAB
#include <mat.h>
#endif

#include <cstdlib>
#include <cassert>
#include <stdlib.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <mpi.h>

#include "include/el/easylogging++.h"
#include <boost/program_options.hpp>

#include "include/spin/Spin.h"
#include "include/spin/SpinCluster.h"
#include "include/spin/SpinClusterAlgorithm.h"
#include "include/spin/SpinCollection.h"
#include "include/spin/SpinData.h"
#include "include/spin/SpinInteraction.h"
#include "include/spin/SpinInteractionComponent.h"
#include "include/spin/SpinInteractionDefine.h"
#include "include/spin/SpinSource.h"
#include "include/spin/SpinState.h"
#include "include/spin/SpinClusterFromLattice.h"

#include "include/kron/KronProd.h"
#include "include/misc/misc.h"
#include "include/misc/xmlreader.h"
#include "include/misc/print_program_options.h"

#include "include/quantum/HilbertSpaceOperator.h"
#include "include/quantum/LiouvilleSpaceOperator.h"
#include "include/quantum/MixedState.h"
#include "include/quantum/PureState.h"
#include "include/quantum/QuantumEvolution.h"
#include "include/quantum/QuantumEvolutionAlgorithm.h"
#include "include/quantum/QuantumOperator.h"
#include "include/quantum/QuantumState.h"

using namespace std;
using namespace arma;

#endif
