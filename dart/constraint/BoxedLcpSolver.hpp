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

#ifndef DART_CONSTRAINT_BOXEDLCPSOLVER_HPP_
#define DART_CONSTRAINT_BOXEDLCPSOLVER_HPP_

#include <Eigen/Core>

namespace dart {
namespace constraint {

/*! BlcpSolver
 *
 *  The boxed linear complementarity problem (BLCP) is defined by
 *
 *  Find \f$(x, y)\f$ such that:
 *    \f{equation*}{
 *    \begin{cases}
 *    A \ x + b = y \\
 *    0 \le x \perp y \ge 0
 *    \end{cases},
 *    \f}
 *  where \f$ x, y, b\f$ are vectors of size \f$n\f$ and \f$ A \f$ is a
 *  \f$n\times n\f$ matrix.
 */
class BoxedLcpSolver
{
public:
  /// Destructor
  virtual ~BoxedLcpSolver() = default;

  /// Solves constriant impulses for a constrained group
  virtual void solve(
      int n,
      double* A,
      double* x,
      double* b,
      int nub,
      double* lo,
      double* hi,
      int* findex) = 0;
  // Note: The function signature is ODE specific for now. Consider changing
  // this to Eigen friendly version once own Dantzig LCP solver is available.

  virtual void solve(
      const Eigen::MatrixXd& A,
      Eigen::VectorXd& x,
      const Eigen::VectorXd& b,
      int nub,
      const Eigen::VectorXd& lo,
      const Eigen::VectorXd& hi) = 0;

  virtual bool canSolve(int n, const double* A) = 0;
};

class LcpErrorMetric
{
public:
  bool canTermninate() const
  {
    return true;
  }
};

} // namespace constraint
} // namespace dart

#endif // DART_CONSTRAINT_BOXEDLCPSOLVER_HPP_

