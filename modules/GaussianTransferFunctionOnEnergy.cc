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


#include <momemta/Logging.h>
#include <momemta/Module.h>
#include <momemta/ParameterSet.h>
#include <momemta/Types.h>
#include <momemta/Math.h>

#include <Math/DistFunc.h>

/** \brief Helper class for Gaussian transfer function modules
 *
 * Base class helping to define TF modules having different behaviours (allowing either to integrate over a TF, or simply evaluate it).
 *
 * \sa GaussianTransferFunctionOnEnergy
 * \sa GaussianTransferFunctionOnEnergyEvaluator
 */
class GaussianTransferFunctionOnEnergyBase: public Module {
    public:

        GaussianTransferFunctionOnEnergyBase(PoolPtr pool, const ParameterSet& parameters): Module(pool, parameters.getModuleName()) {
            m_reco_input = get<LorentzVector>(parameters.get<InputTag>("reco_particle"));

            m_sigma = parameters.get<double>("sigma", 0.10);
            m_sigma_range = parameters.get<double>("sigma_range", 5);
        }

    protected:
        double m_sigma;
        double m_sigma_range;

        // Input
        Value<LorentzVector> m_reco_input;
};

/** \brief Integrate over a transfer function on energy described by a Gaussian distribution
 *
 * This module takes as inputs a LorentzVector and a phase-space point, generates
 * a new LorentzVector with a different energy (keeping direction and invariant mass),
 * and evaluates the transfer function on the "reconstructed" and "generated" energies.
 *
 * The transfer function (TF) is a Gaussian distribution that describes the difference between 
 * the reconstructed and the generated energy (\f$E_{rec}-E_{gen}\f$). The width of the distribution, parametrised as a fraction of \f$E_{gen}\f$, is set as parameter. 
 *
 * The range of the integration is determined using the width of the Gaussian at \f$E_{rec}\f$, integrating over a user-defined 'number of sigmas' `n`: \f$E_{gen} \in \pm n \cdot \sigma \cdot E_{rec}\f$.
 *
 * ### Integration dimension
 *
 * This module requires **1** phase-space point.
 *
 * ### Parameters
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `sigma` | double | Fraction of the energy yielding the width of the Gaussian distribution (with `sigma` at `0.1`, \f$\sigma_{gauss} = 0.1 \cdot E_{gen}\f$). |
 *   | `sigma_range` | double | Range of integration expressed in number of sigma. |
 * 
 * ### Inputs
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `ps_point` | double | Phase-space point generated by CUBA. |
 *   | `reco_particle` | LorentzVector | Input LorentzVector (experimentally reconstructed particle). |
 *
 * ### Outputs
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `output` | LorentzVector | Output *generated* LorentzVector, only differing from *reco_particle* by its energy. |
 *   | `TF_times_jacobian` | double | Product of the TF evaluated on the *reco* and *gen* energies, times the jacobian of the transformation needed stretch the integration range from \f$[0,1]\f$ to the width of the TF, times the jacobian \f$dE/d|P|\f$ due to the fact that the integration is done w.r.t \f$|P|\f$, while the TF is parametrised in terms of energy. |
 * 
 * \ingroup modules
 * \sa GaussianTransferFunctionOnEnergyEvaluator
 */
class GaussianTransferFunctionOnEnergy: public GaussianTransferFunctionOnEnergyBase {
    public:
        GaussianTransferFunctionOnEnergy(PoolPtr pool, const ParameterSet& parameters): GaussianTransferFunctionOnEnergyBase(pool, parameters) {
            m_ps_point = get<double>(parameters.get<InputTag>("ps_point"));
        }

        virtual Status work() override {
            // Estimate the width over which to integrate using the width of the TF at E_rec ...
            const double sigma_E_rec = m_reco_input->E() * m_sigma;

            double range_min = std::max(0., m_reco_input->E() - (m_sigma_range * sigma_E_rec));
            double range_max = m_reco_input->E() + (m_sigma_range * sigma_E_rec);
            double range = (range_max - range_min);

            double gen_E = range_min + range * (*m_ps_point);
            double gen_pt = std::sqrt(SQ(gen_E) - SQ(m_reco_input->M())) / std::cosh(m_reco_input->Eta());

            output->SetCoordinates(
                    gen_pt * std::cos(m_reco_input->Phi()),
                    gen_pt * std::sin(m_reco_input->Phi()),
                    gen_pt * std::sinh(m_reco_input->Eta()),
                    gen_E);

            // ... but compute the width of the TF at E_gen!
            const double sigma_E_gen = gen_E * m_sigma;

            // Compute TF*jacobian, where the jacobian includes the transformation of [0,1]->[range_min,range_max] and d|P|/dE
            *TF_times_jacobian = ROOT::Math::normal_pdf(gen_E, sigma_E_gen, m_reco_input->E()) * range * dP_over_dE(*output);

            return Status::OK;
        }

    private:
        // Input
        Value<double> m_ps_point;

        // Outputs
        std::shared_ptr<LorentzVector> output = produce<LorentzVector>("output");
        std::shared_ptr<double> TF_times_jacobian = produce<double>("TF_times_jacobian");

};

/** \brief Evaluate a transfer function on energy described by a Gaussian distribution
 *
 * This module takes as inputs two LorentzVectors: a 'gen-level' particle (which may be computed using for instance a Block or a 'real' transfer function) and a 'reco-level' particle (experimentally reconstructed). 
 * Assuming the LorentzVectors differ only by their energy, this module returns the value of a transfer function (TF) evaluated on their respective energies.
 *
 * The TF is a Gaussian distribution that describes the difference between the reconstructed and the generated energy (\f$E_{rec}-E_{gen}\f$). The width of the distribution, parametrised as a fraction of \f$E_{gen}\f$, is set as parameter. 
 *
 * ### Integration dimension
 *
 * This module requires **0** phase-space points.
 *
 * ### Parameters
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `sigma` | double | Fraction of the energy yielding the width of the Gaussian distribution (with `sigma` at `0.1`, \f$\sigma_{gauss} = 0.1 \cdot E_{gen}\f$). |
 * 
 * ### Inputs
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `reco_particle` | LorentzVector | Experimentally reconstructed particle. |
 *   | `gen_particle` | LorentzVector | Gen-level particle. |
 *
 * ### Outputs
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `TF` | double | Value of the TF evaluated on the *reco* and *gen* energies. |
 * 
 * \ingroup modules
 * \sa GaussianTransferFunctionOnEnergy
 */
class GaussianTransferFunctionOnEnergyEvaluator: public GaussianTransferFunctionOnEnergyBase {
    public:
        GaussianTransferFunctionOnEnergyEvaluator(PoolPtr pool, const ParameterSet& parameters): GaussianTransferFunctionOnEnergyBase(pool, parameters) {
            m_gen_input = get<LorentzVector>(parameters.get<InputTag>("gen_particle"));
        }

        virtual Status work() override {
            // Compute TF value
            *TF_value = ROOT::Math::normal_pdf(m_gen_input->E(), m_gen_input->E() * m_sigma, m_reco_input->E());

            return Status::OK;
        }

    private:
        // Input
        Value<LorentzVector> m_gen_input;

        // Outputs
        std::shared_ptr<double> TF_value = produce<double>("TF");

};

REGISTER_MODULE(GaussianTransferFunctionOnEnergy);
REGISTER_MODULE(GaussianTransferFunctionOnEnergyEvaluator);
