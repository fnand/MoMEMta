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
#include <TLorentzVector.h>

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
    //
     
    TLorentzVector blor  = TLorentzVector();
    blor.SetPtEtaPhiE(45511.4765625, 2.0693726539611816, 1.9048484563827515, 183206.515625);
    Particle b1 { "bjet1", LorentzVector(blor.Px(), blor.Py(),blor.Pz(), blor.E()), 5 };  
  // Muon
    TLorentzVector mlor  = TLorentzVector();
    mlor.SetPtEtaPhiE(37282.54296875, 0.8650407195091248, -2.555448055267334, 52123.66015625);
    Particle muon { "muon", LorentzVector(mlor.Px(), mlor.Py(),mlor.Pz(), mlor.E()), 13 };
    // Anti b-quark
    //
    TLorentzVector llor  = TLorentzVector();
    llor.SetPtEtaPhiE(30431.296875, 1.264741063117981, -1.1504430770874023, 58653.9765625);
    Particle lj { "ljet1", LorentzVector(llor.Px(), llor.Py(),llor.Pz(), llor.E()), 0 };
    // MET
    TLorentzVector metlor  = TLorentzVector();
    metlor.SetPtEtaPhiE(49928.453125, 0, -0.4170082211494446, 0 );
    LorentzVector met { metlor.Px(), metlor.Py(),metlor.Pz(), metlor.E()};
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
