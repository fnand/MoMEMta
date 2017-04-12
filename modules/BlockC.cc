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


#include <momemta/ParameterSet.h>
#include <momemta/Module.h>
#include <momemta/Solution.h>
#include <momemta/Types.h>
#include <momemta/Math.h>

/** \brief \f$\require{cancel}\f$ Final (main) Block C, describing \f$q_1 q_2 \to \f$X + s_{123} (\to p_3 + s_{12} (\to \cancel{p_1} p_2))\f$
 *
 * This Block addresses the change of variables needed to pass from the standard phase-space
 * parametrisation for \f$p_{1 \dots 3} \times \delta^4\f$ to a parametrisation in terms of the two (squared) masses
 * of the intermediate propagators.
 *
 * The integration is performed over \f$s_{12}, s_{123}\f$ with \f$p_{1 \dots 3}\f$ as input. Per integration point,
 * the LorentzVector of the invisible particle, \f$p_1\f$, is computed based on the following set
 * of equations:
 *
 * - \f$s_{12} = (p_1 + p_2)^2\f$
 * - \f$s_{123} = (p_1 + p_2 + p_3)^2\f$
 * - Conservation of momentum (with \f$\vec{p}_T^{branches}\f$ the total transverse momentum of all branches represented):
 *  - \f$p_{1x} + E_3 \sin\theta_3\cos\phi_3 = - p_{Tx}^{branches}\f$
 *  - \f$p_{1y} + E_3 \sin\theta_3\cos\phi_3 = - p_{Ty}^{branches}\f$
 *
 * Up to four solutions are possible for \f$(E1, \alpha)\f$, where \f$\alpha = 2p_1\dotp_2\f$.
 *
 * ### Integration dimension
 *
 * This module requires **0** phase-space point.
 *
 * ### Global parameters
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `energy` | double | Collision energy. |
 *
 * ### Parameters
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `pT_is_met` | bool, default false | Fix \f$\vec{p}_{T}^{tot} = -\vec{\cancel{E_T}}\f$ or \f$\vec{p}_{T}^{tot} = \sum_{i \in \text{ vis}} \vec{p}_i\f$ |
 *
 * ### Inputs
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `s12` <br/> `s123` <br/>  | double | Squared invariant masses of the propagators. Typically coming from a BreitWignerGenerator or NarrowWidthApproximation module. |
 *   | `p1` ... `p3` | LorentzVector | LorentzVectors of the particles used to reconstruct the event according to the above method. |
 *   | `branches` | vector(LorentzVector) | LorentzVectors of all the other particles in the event, taken into account when computing \f$\vec{p}_{T}^{branches}\f$ (if MET is not used), and checking if the solutions are physical. |
 *   | `met` | LorentzVector, default `met::p4` | LorentzVector of the MET |
 *
 * ### Outputs
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `solutions` | vector(Solution) | Solutions of the change of variable. Each solution embed  the LorentzVectors of the invisible particle (ie. \f$(p_1)\f$) and the associated jacobian. These solutions should be fed as input to the Looper module. |
 *
 * \note This block has NOT been validated and is NOT safe to use.
 *
 * \sa Looper module to loop over the solutions of this Block
 *
 * \ingroup modules
 */

class BlockC: public Module {
    public:

        BlockC(PoolPtr pool, const ParameterSet& parameters): Module(pool, parameters.getModuleName()) {
            sqrt_s = parameters.globalParameters().get<double>("energy");
            pT_is_met = parameters.get<bool>("pT_is_met", false);

            s12  = get<double>(parameters.get<InputTag>("s12"));
            s123 = get<double>(parameters.get<InputTag>("s123"));

	    

            m_particles.push_back(get<LorentzVector>(parameters.get<InputTag>("p2")));
            m_particles.push_back(get<LorentzVector>(parameters.get<InputTag>("p3")));
            
            if (parameters.exists("branches")) {
                auto branches_tags = parameters.get<std::vector<InputTag>>("branches");
                for (auto& t: branches_tags)
                    m_branches.push_back(get<LorentzVector>(t));
            }

            // If the met input is specified, get it, otherwise retrieve default
            // one ("met::p4")
            InputTag met_tag;
            if (parameters.exists("met")) {

                met_tag = parameters.get<InputTag>("met");
            } else {
                met_tag = InputTag({"met", "p4"});
            }

            m_met = get<LorentzVector>(met_tag);

        };

        virtual Status work() override {

            solutions->clear();

            // Don't spend time on unphysical corner of the phase-space
            
            //std::cout << "s12:   " << *s12  << std::endl;
            //std::cout << "s123:  " << *s123 << std::endl;
            //std::cout << "ptr to s12" << s12 << std::endl;
 	    if (*s12 >= *s123 || *s12 >= SQ(sqrt_s) || *s123 >= SQ(sqrt_s) )
                return Status::NEXT;

            const LorentzVector& p2 = *m_particles[0];
            const LorentzVector& p3 = *m_particles[1];
           

 

            // pT will be used to fix the transverse momentum of the reconstructed neutrinos
            // We can either enforce momentum conservation by disregarding the MET, ie:
            //  pT = sum of all the visible particles,
            // Or we can fix it using the MET given as input:
            //  pT = -MET
            // In the latter case, it is the user's job to ensure momentum conservation at
            // the matrix element level (by using the Boost module, for instance).

            LorentzVector pT;
            if (pT_is_met) {
                pT = - *m_met;
            } else {
		//Err, do I need to add that other quark
                pT = p2 + p3;
                for (size_t i = 0; i < m_branches.size(); i++) {
                    pT += *m_branches[i];
                }
            }

            // p1x = alpha1 E1 + beta1 ALPHA + gamma1
            // p1y = ...(2)
            // p1z = ...(3)
	    // E3  = ...(4)            

	    const double p2Sq    = p2.Dot(p2);
	    const double cosphi = std::cos(p3.Phi());
            const double sinphi = std::sin(p3.Phi());
            const double costhe = std::cos(p3.Theta());
            const double sinthe = std::sin(p3.Theta());
	    
	    // Denominator that appears in several of the follwing eq.
	    // No need to compute it multiple times
	    const double denom = 2.*(p2.E() -p2.Pz()*costhe - p2.Px()*cosphi*sinthe - p2.Py()*sinthe*sinphi );


            const double alpha1 = 0.;
            const double beta1  = (cosphi*sinthe)/denom;
            const double gamma1 = -(2*p2.E()*pT.Px() - 2*p2.Pz()*pT.Px()*costhe - *s12*cosphi*sinthe + *s123*cosphi*sinthe -2*p2.Px()*pT.Px()*cosphi*sinthe -2*p2.Py()*pT.Px()*sinthe*sinphi)/denom;


            const double alpha2 = 0.;
            const double beta2  = (sinthe*sinphi)/denom;
            const double gamma2 = -(2*p2.E()*pT.Py() - 2*p2.Pz()*pT.Py()*costhe - 2*p2.Px()*pT.Py()*cosphi*sinthe - 2*p2.Py()*pT.Py()*sinthe*sinphi - *s12*sinthe*sinphi + *s123*sinthe*sinphi)/denom;


            const double alpha3 = 2*(-SQ(p2.E()) + p2.E()*p2.Pz()*costhe + p2.E()*p2.Px()*cosphi*sinthe + p2.E()*p2.Py()*sinthe*sinphi)/(-p2.Pz()*denom);
            const double beta3  = (p2.Px()*cosphi*sinthe + p2.Py()*sinthe*sinphi)/(-p2.Pz()*denom);

            const double gamma3 = (-p2Sq*p2.E() - 2*p2.E()*p2.Px()*pT.Px() - 2*p2.E()*p2.Py()*pT.Py() + p2.E()*(*s12) + p2Sq*p2.Pz()*costhe + 2*p2.Px()*p2.Pz()*pT.Px()*costhe 
				  + 2*p2.Py()*p2.Pz()*pT.Py()*costhe - p2.Pz()*(*s12)*costhe + p2Sq*p2.Px()*cosphi*sinthe   
				  + 2*SQ(p2.Px())*pT.Px()*cosphi*sinthe + 2*p2.Px()*p2.Py()*pT.Py()*cosphi*sinthe - p2.Px()*(*s123)*cosphi*sinthe + p2Sq*p2.Py()*sinthe*sinphi  
                                  + 2*p2.Px()*p2.Py()*pT.Px()*sinthe*sinphi + 2*SQ(p2.Py())*pT.Py()*sinthe*sinphi - p2.Py()*(*s123)*sinthe*sinphi)/(-p2.Pz()*denom);

            const double alpha4 = 0.;
            const double beta4  = -1./denom;
            const double gamma4 = ( *s123- *s12)/denom;



	   


            // a11 E1^2 + a22 ALPHA^2 + a12 E1*ALPHA + a10 E1 + a01 ALPHA + a00 = 0
            // id. with bij

            const double a11 = alpha4*(1- alpha1*sinthe*cosphi - alpha2*sinthe*sinphi - alpha3*costhe);
            const double a22 = beta4*(  -  beta1*sinthe*cosphi -  beta2*sinthe*sinphi - beta3*costhe );
            const double a12 = beta4  - (alpha1*beta4  + alpha4*beta1)*sinthe*cosphi  - (alpha2*beta4  + alpha4*beta2)*sinthe*sinphi  - (alpha3*beta4  + alpha4*beta3)*costhe;
            const double a10 = gamma4 - (alpha1*gamma4 + alpha4*gamma1)*sinthe*cosphi - (alpha2*gamma4 + alpha4*gamma2)*sinthe*sinphi - (alpha3*gamma4 + alpha4*gamma3)*costhe;
            const double a01 = -0.5   - (beta1*gamma4  + beta4*gamma1)*sinthe*cosphi  - (beta2*gamma4  + beta4*gamma2)*sinthe*sinphi  - (beta3*gamma4  + beta4*gamma3)*costhe;
            const double a00 = gamma4* (- gamma1*sinthe*cosphi - gamma2*sinthe*sinphi - gamma3*costhe);

            //std::cout << "a11: " << a11 << "  a22: " << a22 << "  a12: " << a12 << "  a10: " << a10 << "  a01: " << a01 << "  a00: " << a00 << std::endl;
            

            const double b11 = SQ(alpha1) + SQ(alpha2) + SQ(alpha3) - 1;
            const double b22 = SQ(beta1)  + SQ(beta2)  + SQ(beta3);
            const double b12 = 2.*( alpha1*beta1  + alpha2*beta2  + alpha3*beta3  );
            const double b10 = 2.*( alpha1*gamma1 + alpha2*gamma2 + alpha3*gamma3 );
            const double b01 = 2.*( beta1*gamma1  + beta2*gamma2  + beta3*gamma3  );
            const double b00 = SQ(gamma1) + SQ(gamma2) + SQ(gamma3);


            //std::cout << "b11: " << b11 << "  b22: " << b22 << "  b12: " << b12 << "  b10: " << b10 << "  b01: " << b01 << "  b00: " << b00 << std::endl;

            // Find the intersection of the 2 conics (at most 4 real solutions for (E1,ALPHA))
            std::vector<double> E1, ALPHA;
	    //std::cout << "Beg" << std::endl;
            solve2Quads(a11, a22, a12, a10, a01, a00, b11, b22, b12, b10, b01, b00, E1, ALPHA, false);
	    //std::cout << "End" << std::endl;



            // For each solution (E1,ALPHA), find the neutrino 4-momentum p1

            if (E1.size() == 0)
                return Status::NEXT;


            for(unsigned int i=0; i<E1.size(); i++){
                const double e1  = E1.at(i);
                const double alp = ALPHA.at(i);

	
		//FIXME this needs to be changes. What does the alpha requirement translate to?
                if (e1 < 0.) //|| alp < 0.)
                    continue;
		
		// Make sure E3 is larger than zero 
		if (alpha4*e1 + beta4*alp + gamma4  < 0.) //|| alp < 0.)
                    continue;


                LorentzVector p1(
                        alpha1*e1 + beta1*alp + gamma1,
                        alpha2*e1 + beta2*alp + gamma2,
                        alpha3*e1 + beta3*alp + gamma3,
                        e1);

		LorentzVector p3(
                        (alpha4*e1 + beta4*alp + gamma4)*sinthe*cosphi,
                        (alpha4*e1 + beta4*alp + gamma4)*sinthe*sinphi,
                        (alpha4*e1 + beta4*alp + gamma4)*costhe,
                         alpha4*e1 + beta4*alp + gamma4);


	//	std::cout << "Mass of p3 " << p3.M() << std::endl;
		//std::cout << "E3 " << alpha4*e1 + beta4*alp + gamma4 << std::endl;


		//std::cout << "met.Px " << (*m_met).Px() << "  met.Py " << (*m_met).Py() << std:: endl;
		//std::cout << "p1.Px " << p1.Px() << "  p1.Py " << p1.Py() <<  std:: endl;
		//std::cout << "s12       " << *s12 << std::endl;
		//std::cout << "(p1 + p2)^2:  " << p1.Dot(p1) + p2.Dot(p2) + 2*p1.Dot(p2)  << std::endl;
	

		//std::cout << "s123      " << *s123 << std::endl;
                //std::cout << "(p1 + p2 + p3)^2:  " << p1.Dot(p1) + p2.Dot(p2) + p3.Dot(p3) + 2*p1.Dot(p2) + 2*p1.Dot(p3) + 2*p2.Dot(p3)  << std::endl;	
                // Check if solutions are physical
                //
                //
                //std::cout << "p2x + p3x  "  <<p2.Px() +  p3.Px() << std::endl;
		//std::cout << "-pTx  " << - pT.Px()  << std::endl;


                LorentzVector tot = p1 + p2 + p3;
                for (size_t i = 0; i < m_branches.size(); i++) {
                    tot += *m_branches[i];
                }
                double q1Pz = std::abs(tot.Pz() + tot.E()) / 2.;
                double q2Pz = std::abs(tot.Pz() - tot.E()) / 2.;
                if(q1Pz > sqrt_s/2 || q2Pz > sqrt_s/2)
                    continue;

                double jacobian = computeJacobian(p1, p2, p3);
                
		// FIXME do I need to multiply jac by pi*constant?
		Solution s {{p1}, M_PI*jacobian, true };
                solutions->push_back(s);
            }

            return solutions->size() > 0 ? Status::OK : Status::NEXT;
        }

        double computeJacobian(const LorentzVector& p1, const LorentzVector& p2, const LorentzVector& p3) {

            const double E1  = p1.E();
            const double p1x = p1.Px();
            const double p1y = p1.Py();
            const double p1z = p1.Pz();

            const double E2  = p2.E();
            const double p2x = p2.Px();
            const double p2y = p2.Py();
            const double p2z = p2.Pz();

            const double E3  = p3.E();
            //const double p3x = p3.Px();
            //const double p3y = p3.Py();
            //const double p3z = p3.Pz();
	    const double p3phi = p3.Phi();
	    const double p3the = p3.Theta();

                       //
            //
            

	    const double cosphi = std::cos(p3phi);
	    const double sinphi = std::sin(p3phi);
	    const double costhe = std::cos(p3the);
	    const double sinthe = std::sin(p3the);
	    
 
	    const double A  = cosphi*sinthe*(p1x*p2z - p1z*p2x + costhe*(E1*p2x - E2*p1x));
	    const double B  = sinphi*sinthe*(p1z*p2y + p1y*p2z);
	    const double C  = SQ(cosphi)*SQ(sinthe)*(E2*p1z - E1*p2z);              
	    const double D  = SQ(sinphi)*SQ(sinthe)*(E2*p1z - E1*p2z); 
	    const double chi= 2*(p3.Dot(p1) + p3.Dot(p2))/E3;
	    const double F  = SQ(sqrt_s)/(SQ(E3)*E1*sinthe);
	    


       	    double inv_jac = F*std::abs(chi*(E2*p1z - E1*p2z) + 2*E2*(A + B + C + D)); 


            return 1. /inv_jac;
        }

    private:
        double sqrt_s;
        bool pT_is_met;

        // Inputs
        Value<double> s12;
        Value<double> s123;
        std::vector<Value<LorentzVector>> m_particles;
        std::vector<Value<LorentzVector>> m_branches;
        Value<LorentzVector> m_met;

        // Outputs
        std::shared_ptr<SolutionCollection> solutions = produce<SolutionCollection>("solutions");
};
REGISTER_MODULE(BlockC);
