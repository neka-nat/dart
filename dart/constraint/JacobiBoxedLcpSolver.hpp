/*
 * Copyright (c) 2011-2017, The DART development contributors
 * All rights reserved.
 *
 * The list of contributors can be found at:
 *   https://github.com/dartsim/dart/blob/master/LICENSE
 *
 * This file is provided under the following "BSD-style" License:
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef DART_CONSTRAINT_JACOBIBOXEDLCPSOLVER_HPP_
#define DART_CONSTRAINT_JACOBIBOXEDLCPSOLVER_HPP_

#include "dart/constraint/BoxedLcpSolver.hpp"

namespace dart {
namespace constraint {

/// Implementation of projected Gauss-Seidel (PGS) LCP solver.
class JacobiBoxedLcpSolver : public BoxedLcpSolver
{
public:
  struct Option
  {
    int mMaxIteration;
    double sor_w;
    double eps_ea;
    double eps_res;
    double eps_div;

    Option(
        int mMaxIteration = 30,
        double sor_w = 0.9,
        double eps_ea = 10e-3,
        double eps_res = 10e-6,
        double eps_div = 10e-9);
  };

  // Documentation inherited.
  void solve(int n,
      double* A,
      double* x,
      double* b,
      int nub,
      double* lo,
      double* hi,
      int* findex) override;

  // Documentation inherited.
  bool canSolve(int n, double* A) override;

  void setOption(const Option& option);

  const Option& getOption() const;

protected:
  Option mOption;
};

} // namespace constraint
} // namespace dart

#endif // DART_CONSTRAINT_JACOBIBOXEDLCPSOLVER_HPP_

