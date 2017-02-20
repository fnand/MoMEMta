/*
 *  MoMEMta: a modular implementation of the Matrix Element Method
 *  Copyright (C) 2016  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <momemta/ConfigurationReader.h>
#include <momemta/Logging.h>
#include <momemta/MoMEMta.h>
#include <momemta/Unused.h>

#include <TH1D.h>

#include <chrono>

using namespace std::chrono;
using namespace momemta;

int main(int argc, char** argv) {

    UNUSED(argc);
    UNUSED(argv);

    logging::set_level(logging::level::debug);

    ParameterSet lua_parameters;
    lua_parameters.set("USE_TF", true);
    lua_parameters.set("USE_PERM", true);

    ConfigurationReader configuration("../examples/VLQ_test.lua", lua_parameters);

    // Change top mass
    configuration.getGlobalParameters().set("top_mass", 173.);

    MoMEMta weight(configuration.freeze());

    // Electron
    //Particle electron { "electron", LorentzVector(16.171895980835, -13.7919054031372, -3.42997527122497, 21.5293197631836), -11 };
    // b-quark
//    Particle b1 { "bjet1", LorentzVector(-26.7908325195313, -30.59294128418, 140.144721984863, 146.66259765625), 5 };
    Particle b1 { "bjet1", LorentzVector(-26790.8325195313, -30592.94128418, 140144.721984863, 146662.59765625), 5 };  
  // Muon
    //Particle muon { "muon", LorentzVector(35., 17.0896110534668, 53., 66.), -13 };
    Particle muon { "muon", LorentzVector(35000., 17089.6110534668, 53000., 66000.), -13 };
    // Anti b-quark
    //Particle lj { "ljet1", LorentzVector(26.6368963276919, -68.04161840859496, 718.2086419584481, 721.955125), 1 };
    Particle lj { "ljet1", LorentzVector(26636.8963276919, -68041.61840859496, 718208.6419584481, 721955.125), 1 };
    // MET
    //LorentzVector met { -40.60706386098534, 79.58668685644121 , 0 , 0};
    LorentzVector met { -40607.06386098534, 79586.68685644121 , 0 , 0};
    auto start_time = system_clock::now();
    std::vector<std::pair<double, double>> weights = weight.computeWeights({muon, b1, lj}, met);
    auto end_time = system_clock::now();

    LOG(debug) << "Result:";
    for (const auto& r: weights) {
        LOG(debug) << r.first << " +- " << r.second;
    }

    LOG(debug) << "Integration status: " << (int) weight.getIntegrationStatus();

    InputTag dmemInputTag {"dmem", "hist"};
    bool exists = weight.getPool().exists(dmemInputTag);

    LOG(debug) << "Hist in pool: " << exists;

    if (exists) {
        Value<TH1D> dmem = weight.getPool().get<TH1D>(dmemInputTag);
        LOG(debug) << "DMEM integral: " << dmem->Integral();
    }

    LOG(info) << "Weight computed in " << std::chrono::duration_cast<milliseconds>(end_time - start_time).count() << "ms";


    return 0;
}
