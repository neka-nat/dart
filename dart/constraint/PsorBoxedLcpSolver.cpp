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

#include "dart/constraint/PsorBoxedLcpSolver.hpp"

#include <cmath>
#include "dart/external/odelcpsolver/lcp.h"
#include "dart/math/Constants.hpp"

#define PGS_EPSILON 10e-9

namespace dart {
namespace constraint {

//==============================================================================
PsorBoxedLcpSolver::Option::Option(
    int itermax, double sor_w, double eps_ea, double eps_res, double eps_div)
  : mMaxIteration(itermax),
    sor_w(sor_w),
    eps_ea(eps_ea),
    eps_res(eps_res),
    eps_div(eps_div)
{
  // Do nothing
}

//==============================================================================
void PsorBoxedLcpSolver::solve(int n,
    double* A,
    double* x,
    double* b,
    int /*nub*/,
    double* lo,
    double* hi,
    int* findex)
{
  const int nskip = dPAD(n);

  // LDLT solver will work !!!
  //if (nub == n)
  //{
  //	return LDLTSolver(n,nskip,A,x,b)
  //}

  int i, j, iter, idx, n_new;
  bool sentinel;
  double old_x, new_x, hi_tmp, lo_tmp, dummy, ea;
  double * A_ptr;
  double one_minus_sor_w = 1.0 - (mOption.sor_w);

  //--- ORDERING & SCALING & INITIAL LOOP & Test
  int* order = new int[n];

  n_new = 0;
  sentinel = true;
  for (i = 0 ; i < n ; i++)
  {
    // ORDERING
    if ( A[nskip*i + i] < mOption.eps_div )
    {
      x[i] = 0.0;
      continue;
    }
    order[n_new++] = i;

    // INITIAL LOOP
    A_ptr = A + nskip*i;
    new_x = b[i];
    old_x = x[i];

    for (j = 0 ; j < i ; j++)
      new_x -= A_ptr[j]*x[j];
    for (j = i + 1 ; j < n ; j++)
      new_x -= A_ptr[j]*x[j];

    new_x = new_x/A[nskip*i + i];

    if (findex[i] >= 0)	// friction index
    {
      hi_tmp = hi[i] * x[findex[i]];
      lo_tmp = -hi_tmp;

      if (new_x > hi_tmp)
        x[i] = hi_tmp;
      else if (new_x < lo_tmp)
        x[i] = lo_tmp;
      else
        x[i] = new_x;
    }
    else					// no friction index
    {
      if (new_x > hi[i])
        x[i] = hi[i];
      else if (new_x < lo[i])
        x[i] = lo[i];
      else
        x[i] = new_x;
    }

    // TEST
    if (sentinel)
    {
      ea = std::abs(x[i] - old_x);
      if (ea > mOption.eps_res)
        sentinel = false;
    }
  }
  if (sentinel)
  {
    delete[] order;
    // return true;
    return;
  }

  // SCALING
  for (i = 0 ; i < n_new ; i++)
  {
    idx = order[i];

    dummy = 1.0/A[nskip*idx + idx];  // diagonal element
    b[idx] *= dummy;
    for (j = 0 ; j < n ; j++)
      A[nskip*idx + j] *= dummy;
  }

  //--- ITERATION LOOP
  for (iter = 1 ; iter < mOption.mMaxIteration; iter++)
  {
    //--- RANDOMLY_REORDER_CONSTRAINTS
#if LCP_PGS_RANDOMLY_REORDER_CONSTRAINTS
    if ((iter & 7)==0)
    {
      int tmp, swapi;
      for (i = 1 ; i < n_new ; i++)
      {
        tmp = order[i];
        swapi = dRandInt(i+1);
        order[i] = order[swapi];
        order[swapi] = tmp;
      }
    }
#endif

    sentinel = true;

    //-- ONE LOOP
    for (i = 0 ; i < n_new ; i++)
    {
      idx = order[i];

      A_ptr = A + nskip*idx;
      new_x = b[idx];
      old_x = x[idx];

      for (j = 0 ; j < idx ; j++)
        new_x -= A_ptr[j]*x[j];
      for (j = idx + 1 ; j < n ; j++)
        new_x -= A_ptr[j]*x[j];

      new_x = (mOption.sor_w * new_x) + (one_minus_sor_w * old_x);

      if (findex[idx] >= 0)	// friction index
      {
        hi_tmp = hi[idx] * x[findex[idx]];
        lo_tmp = -hi_tmp;

        if (new_x > hi_tmp)
          x[idx] = hi_tmp;
        else if (new_x < lo_tmp)
          x[idx] = lo_tmp;
        else
          x[idx] = new_x;
      }
      else					// no friction index
      {
        if (new_x > hi[idx])
          x[idx] = hi[idx];
        else if (new_x < lo[idx])
          x[idx] = lo[idx];
        else
          x[idx] = new_x;
      }

      if ( sentinel && std::abs(x[idx]) > mOption.eps_div)
      {
        ea = std::abs((x[idx] - old_x)/x[idx]);
        if (ea > mOption.eps_ea)
          sentinel = false;
      }
    }

    if (sentinel)
      break;
  }
  delete[] order;

  //return sentinel; // TODO(JS): Consider providing the result passing sentinel
  // in it.
}

//==============================================================================
void PsorBoxedLcpSolver::solve(const Eigen::MatrixXd& A,
    Eigen::VectorXd& x,
    const Eigen::VectorXd& b,
    int nub,
    const Eigen::VectorXd& lo,
    const Eigen::VectorXd& hi)
{
  A;
  x;
  b;
  nub;
  lo;
  hi;
}

//==============================================================================
bool PsorBoxedLcpSolver::canSolve(int n, const double* A)
{
  const int nskip = dPAD(n);

  // Return false if A has zero-diagonal or A is nonsymmetric matrix
  for (auto i = 0; i < n; ++i)
  {
    if (A[nskip * i + i] < PGS_EPSILON)
      return false;

    for (auto j = 0; j < n; ++j)
    {
      if (std::abs(A[nskip * i + j] - A[nskip * j + i]) > PGS_EPSILON)
        return false;
    }
  }

  return true;
}

//==============================================================================
void PsorBoxedLcpSolver::setOption(const PsorBoxedLcpSolver::Option& option)
{
  mOption = option;
}

//==============================================================================
const PsorBoxedLcpSolver::Option& PsorBoxedLcpSolver::getOption() const
{
  return mOption;
}

}  // namespace constraint
}  // namespace dart
