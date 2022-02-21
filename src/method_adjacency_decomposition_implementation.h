
//-------------------------------------------------------------------------------//
// Author: Stefan Lörwald, Universität Heidelberg                                //
// License: CC BY-NC 4.0 http://creativecommons.org/licenses/by-nc/4.0/legalcode //
//-------------------------------------------------------------------------------//

#pragma once

#include <tuple>

#include "maps.h"
#include "matrix.h"
#include "names.h"
#include "tags.h"

namespace panda
{
   namespace implementation
   {
      /// Helper function for adjacency decomposition.
      template <template <typename, typename> class JobManagerType, typename Integer, typename TagType>
      void adjacencyDecomposition(int, char**, const std::tuple<Matrix<Integer>, Names, Maps, Matrix<Integer>>&, TagType);

       template<typename Integer, typename TagType>
       void print_class_representatives(const Matrix<Integer> &rows, const Maps &maps, TagType);
   }
}

#include "method_adjacency_decomposition_implementation.eti"

