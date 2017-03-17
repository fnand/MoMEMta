// 
// *  This file was automatically generated by MoMEMta-MaGMEE,
// *  A MadGraph Matrix Element Exporter plugin for MoMEMta.
// *
// *  It is subject to MoMEMta-MaGMEE's license and copyright:
// *
// *  Copyright (C) 2016  Universite catholique de Louvain (UCL), Belgium
// *
// *  This program is free software: you can redistribute it and/or modify
// *  it under the terms of the GNU General Public License as published by
// *  the Free Software Foundation, either version 3 of the License, or
// *  (at your option) any later version.
// *
// *  This program is distributed in the hope that it will be useful,
// *  but WITHOUT ANY WARRANTY; without even the implied warranty of
// *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// *  GNU General Public License for more details.
// *
// *  You should have received a copy of the GNU General Public License
// *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
// 

#include <string> 
#include <utility> 
#include <vector> 
#include <map> 

#include <P1_Sigma_sm_no_b_mass_dux_bvtxtambx.h> 
#include <HelAmps_sm_no_b_mass.h> 

#include <momemta/ParameterSet.h> 
#include <momemta/SLHAReader.h> 

namespace SgTop_sm_no_b_mass 
{

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: d u~ > b t~ WEIGHTED<=4 / h z @1
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > vt~ ta- WEIGHTED<=2
// Process: s c~ > b t~ WEIGHTED<=4 / h z @1
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > vt~ ta- WEIGHTED<=2

//--------------------------------------------------------------------------

// Initialize process.

P1_Sigma_sm_no_b_mass_dux_bvtxtambx::P1_Sigma_sm_no_b_mass_dux_bvtxtambx(const
    ParameterSet& configuration)
{

  std::string param_card = configuration.get < std::string > ("card"); 
  params.reset(new Parameters_sm_no_b_mass(SLHA::Reader(param_card))); 

  // Set external particle masses for this matrix element
  mME.push_back(std::ref(params->ZERO)); 
  mME.push_back(std::ref(params->ZERO)); 
  mME.push_back(std::ref(params->ZERO)); 
  mME.push_back(std::ref(params->ZERO)); 
  mME.push_back(std::ref(params->mdl_MTA)); 
  mME.push_back(std::ref(params->ZERO)); 

  mapFinalStates[{5, -16, 15, -5}] = 
  {
    {
      &P1_Sigma_sm_no_b_mass_dux_bvtxtambx::matrix_1_dux_btx_no_hz_tx_wmbx_wm_vtxtam, 
      true, 
      {
        std::make_pair(1, -2), std::make_pair(3, -4)
      }, 
      64, 
      36
    }
  }; 

}

//--------------------------------------------------------------------------
// Evaluate |M|^2, return a map of final states

std::map < std::pair < int, int > , double >
    P1_Sigma_sm_no_b_mass_dux_bvtxtambx::compute(const std::pair <
    std::vector<double> , std::vector<double> > &initialMomenta, const
    std::vector < std::pair < int, std::vector<double> > > &finalState)
{

  // Set initial particle momenta
  momenta[0] = (double * ) (&initialMomenta.first[0]); 
  momenta[1] = (double * ) (&initialMomenta.second[0]); 

  // Suppose final particles are passed in the "correct" order
  std::vector<int> selectedFinalState(6 - 2); 
  for (size_t index = 0; index < (6 - 2); index++ )
  {
    selectedFinalState[index] = finalState[index].first; 
    momenta[index + 2] = (double * ) (&finalState[index].second[0]); 
  }

  // Set the event specific parameters
  params->updateParameters(); 
  params->updateCouplings(); 

  // Initialise result object
  std::map < std::pair < int, int > , double > result; 

  // Define permutation
  int perm[6]; 
  for(int i = 0; i < 6; i++ )
  {
    perm[i] = i; 
  }

  for(auto &me: mapFinalStates[selectedFinalState])
  {

    double me_sum = 0; 
    double me_mirror_sum = 0; 

    for(int ihel = 0; ihel < 64; ihel++ )
    {

      if(me.goodHel[ihel])
      {

        double sum = 0.; 
        calculate_wavefunctions(perm, helicities[ihel]); 
        double meTemp = me.callback( * this); 
        sum += meTemp; 
        me_sum += meTemp/me.denominator; 

        if(me.hasMirrorProcess)
        {
          perm[0] = 1; 
          perm[1] = 0; 
          // Calculate wavefunctions
          calculate_wavefunctions(perm, helicities[ihel]); 
          // Mirror back
          perm[0] = 0; 
          perm[1] = 1; 
          meTemp = me.callback( * this); 
          sum += meTemp; 
          me_mirror_sum += meTemp/me.denominator; 
        }

        if( !sum)
          me.goodHel[ihel] = false; 
      }
    }

    for (auto const &initialState: me.initialStates)
    {
      result[initialState] = me_sum; 
      if (me.hasMirrorProcess)
        result[std::make_pair(initialState.second, initialState.first)] =
            me_mirror_sum;
    }
  }


  return result; 
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void P1_Sigma_sm_no_b_mass_dux_bvtxtambx::calculate_wavefunctions(const int
    perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  static std::complex<double> w[9][18]; 

  // Calculate all wavefunctions
  ixxxxx(&momenta[perm[0]][0], mME[0], hel[0], +1, w[0]); 
  oxxxxx(&momenta[perm[1]][0], mME[1], hel[1], -1, w[1]); 
  oxxxxx(&momenta[perm[2]][0], mME[2], hel[2], +1, w[2]); 
  ixxxxx(&momenta[perm[3]][0], mME[3], hel[3], -1, w[3]); 
  oxxxxx(&momenta[perm[4]][0], mME[4], hel[4], +1, w[4]); 
  FFV2_3(w[3], w[4], params->GC_100, params->mdl_MW, params->mdl_WW, w[5]); 
  ixxxxx(&momenta[perm[5]][0], mME[5], hel[5], -1, w[6]); 
  FFV2_2(w[6], w[5], params->GC_100, params->mdl_MT, params->mdl_WT, w[7]); 
  FFV2_3(w[0], w[1], params->GC_100, params->mdl_MW, params->mdl_WW, w[8]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV2_0(w[7], w[2], w[8], params->GC_100, amp[0]); 

}
double P1_Sigma_sm_no_b_mass_dux_bvtxtambx::matrix_1_dux_btx_no_hz_tx_wmbx_wm_vtxtam() 
{

  static std::complex<double> ztemp; 
  static std::complex<double> jamp[1]; 
  // The color matrix
  static const double denom[1] = {1}; 
  static const double cf[1][1] = {{9}}; 

  // Calculate color flows
  jamp[0] = -amp[0]; 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(int i = 0; i < 1; i++ )
  {
    ztemp = 0.; 
    for(int j = 0; j < 1; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  return matrix; 
}



}

// Register matrix element with MoMEMta
#include <momemta/MatrixElementFactory.h> 
REGISTER_MATRIX_ELEMENT("SgTop_sm_no_b_mass_P1_Sigma_sm_no_b_mass_dux_bvtxtambx", SgTop_sm_no_b_mass::P1_Sigma_sm_no_b_mass_dux_bvtxtambx); 
