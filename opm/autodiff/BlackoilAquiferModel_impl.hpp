namespace Opm {


    template<typename TypeTag>
    BlackoilAquiferModel<TypeTag>::
    BlackoilAquiferModel(Simulator& ebosSimulator,
                      const ModelParameters& param,
                      const bool terminal_output)
        : ebosSimulator_(ebosSimulator)
        , param_(param)
        , terminal_output_(terminal_output)
        , has_solvent_(GET_PROP_VALUE(TypeTag, EnableSolvent))
        , has_polymer_(GET_PROP_VALUE(TypeTag, EnablePolymer))
    {
        const auto& eclState = ebosSimulator_.gridManager().eclState();
        phase_usage_ = phaseUsageFromDeck(eclState);

        active_.resize(phase_usage_.MaxNumPhases, false);
        for (int p = 0; p < phase_usage_.MaxNumPhases; ++p) {
            active_[ p ] = phase_usage_.phase_used[ p ] != 0;
        }

        const auto& gridView = ebosSimulator_.gridView();

        // calculate the number of elements of the compressed sequential grid. this needs
        // to be done in two steps because the dune communicator expects a reference as
        // argument for sum()
        number_of_cells_ = gridView.size(/*codim=*/0);
        global_nc_ = gridView.comm().sum(number_of_cells_);
        gravity_ = ebosSimulator_.problem().gravity()[2];
        aquifers_ = hack_init(ebosSimulator_);
    }



    // called at the beginning of a time step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: beginTimeStep() 
    {
        // Right now it doesn't do shit.
    }

    // called at the end of a time step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: timeStepSucceeded()
    {
        // Right now it doesn't do shit.
    }

    // called at the beginning of a report step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: beginReportStep(const int time_step) 
    {
        // Right now it doesn't do shit.
    }

    // called at the end of a report step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: endReportStep() 
    {
        // Right now it just spits out the constants for each aquifers
        // We are using the simple integer indexing for the aquifers
        for (int i = 0; i < numAquifers(); ++i)
        {
            std::cout << "Aquifer[" << i << "] -> id " << aquifers()[i].cell_id()
                      << " : Tc = " << aquifers()[i].time_constant()
                      << ", beta = " << aquifers()[i].aquifer_influx_constant() << std::endl;
        }
    }

    // Get the last report step
    template<typename TypeTag>
    const SimulatorReport& 
    BlackoilAquiferModel<TypeTag>:: lastReport() const 
    {
        for (auto i = aquifers_.begin(); i != aquifers_.end(); ++i){
            (*i).print_private_members();
        }
        return last_report_;
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::
    assemble( const SimulatorTimerInterface& timer,
              const int iterationIdx                )
    {
        last_report_ = SimulatorReport();

        if ( ! aquifersActive() ) {
            return;
        }

        
        // We need to update the reservoir pressures connected to the aquifer
        updateConnectionIntensiveQuantities();

        if (iterationIdx == 0) {
            // We can do the Table check and coefficients update in this function
            // For now, it does nothing!
            prepareTimeStep();
        }

        if (iterationIdx == 0) {
            calculateExplicitQuantities();
        }

        if (param_.solve_aquifereq_initially_ && iterationIdx == 0) {
            // solve the aquifer equations as a pre-processing step
            last_report_ = solveAquiferEq(timer);
        }

        assembleAquiferEq(timer);

        last_report_.converged = true;
    }

    // Protected function: Update the primary variables
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: updatePrimaryVariables() 
    {
        // Right now it doesn't do shit.
    }

    // Protected function: Init the primary variables
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: initPrimaryVariablesEvaluation() const 
    {
        // Right now it doesn't do shit.
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: updateConnectionIntensiveQuantities() const
    {
        ElementContext elemCtx(ebosSimulator_);
        const auto& gridView = ebosSimulator_.gridView();
        const auto& elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();
        for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
             elemIt != elemEndIt;
             ++elemIt)
        {
            elemCtx.updatePrimaryStencil(*elemIt);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
        }
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: calculateExplicitQuantities() const
    {
        // For now we can only do one aquifer, so we don't loop over all...
        // for (int aqid; aqid < numAquifers(); aqid++)
        // {
        //     (aquifers.at(aqid))->calculateExplicitQuantities(ebosSimulator_);
        // }

        /*
         for (auto& aquifer : aquifer_container_) {
             aquifer->calculateExplicitQuantities(ebosSimulator_);
         }*/
    }


    template<typename TypeTag>
    SimulatorReport 
    BlackoilAquiferModel<TypeTag>:: solveAquiferEq(const SimulatorTimerInterface& timer)
    {
        // We need to solve the equilibrium equation first to
        // obtain the initial pressure of water in the aquifer
        SimulatorReport report = SimulatorReport();
        return report;
    }

    // Protected function: Return number of components in the model.
    template<typename TypeTag>
    int
    BlackoilAquiferModel<TypeTag>:: numComponents() const
    {
        if (numPhases() == 2) {
            return 2;
        }
        int numComp = FluidSystem::numComponents;
        if (has_solvent_) {
            numComp ++;
        }

        return numComp;
    }

    // Protected function: Return number of aquifers in the model.
    template<typename TypeTag>
    int
    BlackoilAquiferModel<TypeTag>:: numAquifers() const
    {
        return aquifers_.size();
    }

    // Protected function: Return number of phases in the model.
    template<typename TypeTag>
    int
    BlackoilAquiferModel<TypeTag>:: numPhases() const
    {
        // Not implemented yet!!!!!!!!!!!!
        return -1;
        // return aquifers()->number_of_phases;
    }


    // Protected function: returns the phase index in ebos
    template<typename TypeTag>
    int
    BlackoilAquiferModel<TypeTag>:: flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
    {
        const auto& pu = phase_usage_;
        if (active_[Water] && pu.phase_pos[Water] == phaseIdx)
            return FluidSystem::waterPhaseIdx;
        if (active_[Oil] && pu.phase_pos[Oil] == phaseIdx)
            return FluidSystem::oilPhaseIdx;
        if (active_[Gas] && pu.phase_pos[Gas] == phaseIdx)
            return FluidSystem::gasPhaseIdx;

        assert(phaseIdx < 3);
        // for other phases return the index
        return phaseIdx;
    }

    // Protected function which calls the individual aquifer models
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: assembleAquiferEq(const SimulatorTimerInterface& timer)
    {
        // Still doesn't do shit
    }

    // Protected function
    // some preparation work, mostly related to group control and RESV,
    // at the beginning of each time step (Not report step)
    template<typename TypeTag>
    void BlackoilAquiferModel<TypeTag>:: prepareTimeStep()
    {
        // Not implemented yet!
    }

    // Protected function: Returns a reference to the aquifers members in the model
    template<typename TypeTag>
    const std::vector< AquiferCarterTracy<TypeTag> >&
    BlackoilAquiferModel<TypeTag>:: aquifers()
    {
        return aquifers_;
    }

    // Protected function
    /// return true if aquifers are available in the reservoir
    template<typename TypeTag>
    bool
    BlackoilAquiferModel<TypeTag>:: aquifersActive() const
    {
        return aquifers_active_;
    }

    // Protected function
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: setAquifersActive(const bool aquifers_active)
    {
        aquifers_active_ = aquifers_active;
    }

    // Protected function
    /// return true if Aquifers are available on this process
    template<typename TypeTag>
    bool
    BlackoilAquiferModel<TypeTag>:: localAquifersActive() const
    {
        return numAquifers() > 0;
    }

    // Begin the hack to initialize the aquifers in the deck
    template<typename TypeTag>
    std::vector< AquiferCarterTracy<TypeTag> >
    BlackoilAquiferModel<TypeTag>:: hack_init(const Simulator& ebosSimulator)//, std::vector< AquiferCarterTracy<TypeTag> >& aquifers)
    {
        std::vector< AquiferCarterTracy<TypeTag> > aquifers;
        /** Begin hack!!!!! */
        const auto& deck = ebosSimulator.gridManager().deck();
        const auto& eclState = ebosSimulator.gridManager().eclState();
        
        if (!deck.hasKeyword("AQUCT")){
            std::cout << "Nothing is shown! Where is AQUCT!????" << std::endl;
        }

        const auto& aquctKeyword = deck.getKeyword("AQUCT");
        std::vector<AQUCT_params> aquctParams;
        // Resize the parameter vector container based on row entries in aquct
        // We do the same for aquifers too because number of aquifers is assumed to be for each entry in aquct
        aquctParams.resize(aquctKeyword.size());
        // aquifers.resize(aquctKeyword.size());

        const int tableID = 0;

        std::cout << "Parsing AQUCT stuff" << std::endl;
        for (size_t aquctRecordIdx = 0; aquctRecordIdx < aquctKeyword.size(); ++ aquctRecordIdx) 
        {
            const auto& aquctRecord = aquctKeyword.getRecord(aquctRecordIdx);

            aquctParams.at(aquctRecordIdx).aquiferID = aquctRecord.getItem("AQUIFER_ID").template get<int>(0);
            aquctParams.at(aquctRecordIdx).h = aquctRecord.getItem("THICKNESS_AQ").template get<double>(0);
            aquctParams.at(aquctRecordIdx).phi_aq = aquctRecord.getItem("PORO_AQ").template get<double>(0);
            aquctParams.at(aquctRecordIdx).d0 = aquctRecord.getItem("DAT_DEPTH").getSIDouble(0);
            aquctParams.at(aquctRecordIdx).C_t = aquctRecord.getItem("C_T").template get<double>(0);
            aquctParams.at(aquctRecordIdx).r_o = aquctRecord.getItem("RAD").getSIDouble(0);
            aquctParams.at(aquctRecordIdx).k_a = aquctRecord.getItem("PERM_AQ").getSIDouble(0);
            aquctParams.at(aquctRecordIdx).theta = aquctRecord.getItem("INFLUENCE_ANGLE").template get<double>(0);
            aquctParams.at(aquctRecordIdx).c1 = 0.008527; // We are using SI
            aquctParams.at(aquctRecordIdx).c2 = 6.283;
            aquctParams.at(aquctRecordIdx).inftableID = aquctRecord.getItem("TABLE_NUM_INFLUENCE_FN").template get<int>(0) - 1;
            aquctParams.at(aquctRecordIdx).pvttableID = aquctRecord.getItem("TABLE_NUM_WATER_PRESS").template get<int>(0) - 1;

            std::cout << aquctParams.at(aquctRecordIdx).inftableID << std::endl;
            // Get the correct influence table values
            const auto& aqutabTable = eclState.getTableManager().getAqutabTables().getTable(aquctParams.at(aquctRecordIdx).inftableID);
            const auto& aqutab_tdColumn = aqutabTable.getColumn(0);
            const auto& aqutab_piColumn = aqutabTable.getColumn(1);
            aquctParams.at(aquctRecordIdx).td = aqutab_tdColumn.vectorCopy();
            aquctParams.at(aquctRecordIdx).pi = aqutab_piColumn.vectorCopy();

            // We determine the cell perforation here.
            int cellID = 10 + aquctRecordIdx;

            aquctParams.at(aquctRecordIdx).cell_id = cellID;

            // We do not have mu_w as it has to be calculated from pvttable
            aquifers.push_back(Aquifer_object( aquctParams.at(aquctRecordIdx) ));
        }

        // I want to deliberately add another aquifer
        aquifers.push_back( Aquifer_object(99) );

        // aquifers_ = aquifers;

        return aquifers;
    }

} // namespace Opm