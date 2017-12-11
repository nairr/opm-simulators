/*
  Copyright 2017 TNO - Heat Transfer & Fluid Dynamics, Modelling & Optimization of the Subsurface
  Copyright 2017 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_AQUIFERCT_HEADER_INCLUDED
#define OPM_AQUIFERCT_HEADER_INCLUDED


#include <opm/autodiff/AQUCT_params.hpp>
#include <opm/autodiff/BlackoilAquiferModel.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/core/props/BlackoilPhases.hpp>


#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <string>
#include <memory>
#include <vector>
#include <cassert>

namespace Opm
{

    template<typename TypeTag>
    class AquiferCarterTracy
    {

        public:
            typedef BlackoilModelParameters ModelParameters;

            static const int Water = BlackoilPhases::Aqua;
            static const int Oil = BlackoilPhases::Liquid;
            static const int Gas = BlackoilPhases::Vapour;

            typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
            typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
            typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
            typedef typename GET_PROP_TYPE(TypeTag, Indices) BlackoilIndices;
            typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
            typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

            static const int numEq = BlackoilIndices::numEq;
            typedef double Scalar;

            typedef DenseAd::Evaluation<double, /*size=*/numEq> Eval;

            typedef Ewoms::BlackOilPolymerModule<TypeTag> PolymerModule;

            static const bool has_solvent = GET_PROP_VALUE(TypeTag, EnableSolvent);
            static const bool has_polymer = GET_PROP_VALUE(TypeTag, EnablePolymer);
            static const int contiSolventEqIdx = BlackoilIndices::contiSolventEqIdx;
            static const int contiPolymerEqIdx = BlackoilIndices::contiPolymerEqIdx;


            AquiferCarterTracy(int cell_id)
            : phi_aq_ (1.0), //
              C_t_ (1.0), //
              r_o_ (1.0), //
              k_a_ (1.0), //
              c1_ (1.0),
              h_ (1.0), //
              theta_ (1.0), //
              c2_ (1.0), //
              d0_ (1.0),
              cell_idx_ (cell_id)
            {
                mu_w_ = 1000.0;
                aqutab_td_.push_back(1.0);
                aqutab_pi_.push_back(1.0);
                aquiferID_ = 1;
                inftableID_ = 1;
                pvttableID_ = 1;
            }

            explicit AquiferCarterTracy(const AQUCT_params& params)
            : phi_aq_ (params.phi_aq), //
              C_t_ (params.C_t), //
              r_o_ (params.r_o), //
              k_a_ (params.k_a), //
              c1_ (params.c1),
              h_ (params.h), //
              theta_ (params.theta), //
              c2_ (params.c2), //
              d0_ (params.d0),
              aqutab_td_ (params.td),
              aqutab_pi_ (params.pi),
              aquiferID_ (params.aquiferID),
              inftableID_ (params.inftableID),
              pvttableID_ (params.pvttableID),
              cell_idx_ (params.cell_id)
            {
                mu_w_ = 1000.0;
            }

            inline const PhaseUsage&
            phaseUsage() const
            {
                assert(phase_usage_);

                return *phase_usage_;
            }

            inline int
            flowPhaseToEbosCompIdx( const int phaseIdx ) const
            {
                const auto& pu = phaseUsage();
                if (active()[Water] && pu.phase_pos[Water] == phaseIdx)
                    return BlackoilIndices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                if (active()[Oil] && pu.phase_pos[Oil] == phaseIdx)
                    return BlackoilIndices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                if (active()[Gas] && pu.phase_pos[Gas] == phaseIdx)
                    return BlackoilIndices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);

                // for other phases return the index
                return phaseIdx;
            }

            inline int
            flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
            {
                const auto& pu = phaseUsage();
                if (active()[Water] && pu.phase_pos[Water] == phaseIdx) {
                    return FluidSystem::waterPhaseIdx;
                }
                if (active()[Oil] && pu.phase_pos[Oil] == phaseIdx) {
                    return FluidSystem::oilPhaseIdx;
                }
                if (active()[Gas] && pu.phase_pos[Gas] == phaseIdx) {
                    return FluidSystem::gasPhaseIdx;
                }

                assert(phaseIdx < 3);
                // for other phases return the index
                return phaseIdx;
            }

            /* Made into public for testing only!!!!!!. Must be protected */
            inline const Scalar time_constant() const
            {
                Scalar Tc = mu_w_*phi_aq_*C_t_*r_o_*r_o_/(k_a_*c1_);
                return Tc;
            }

            /* Made into public for testing only!!!!!!. Must be protected */
            inline const Scalar aquifer_influx_constant() const
            {
                Scalar beta = c2_*h_*theta_*phi_aq_*C_t_*r_o_*r_o_;
                return beta;
            }

            inline void print_private_members() const
            {
                std::cout << "Aquifer CT #" << aquiferID_ << std::endl;
                auto ita = aqutab_td_.cbegin();
                auto f_lambda = [&ita] (double i) {std::cout << *ita++ << "    " << i << std::endl;};
                std::for_each( aqutab_pi_.cbegin(), aqutab_pi_.cend(), f_lambda );
            }

            /* Made into public for testing only!!!!!!. Must be protected */
            inline const int cell_id() const
            {
                return cell_idx_;
            }



        protected:
            const PhaseUsage* phase_usage_;
            const std::vector<bool>* active_;

            inline const std::vector<bool>& 
            active() const
            {
                assert(active_);
                return *active_;
            }

            // Aquifer ID, and other IDs
            int aquiferID_, inftableID_, pvttableID_;

            // Grid variables
            int cell_idx_;

            // Variables constants
            Scalar mu_w_ , //water viscosity
                   phi_aq_ , //aquifer porosity
                   d0_,     // aquifer datum depth
                   C_t_ , //total compressibility
                   r_o_ , //aquifer inner radius
                   k_a_ , //aquifer permeability
                   c1_, // 0.008527 (METRIC, PVT-M); 0.006328 (FIELD); 3.6 (LAB)
                   h_ , //aquifer thickness
                   theta_ , //angle subtended by the aquifer boundary
                   c2_ ; //6.283 (METRIC, PVT-M); 1.1191 (FIELD); 6.283 (LAB).

            // Variables for influence table
            std::vector<Scalar> aqutab_td_, aqutab_pi_;
            

    }; // class AquiferCarterTracy


} // namespace Opm

#endif