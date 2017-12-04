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
    }

    // called at the beginning of a time step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::
    beginTimeStep() {
        // Right now it doesn't do shit.
    }

    // called at the end of a time step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::
    timeStepSucceeded() {
        // Right now it doesn't do shit.
    }

    // called at the beginning of a report step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::
    beginReportStep(const int time_step) {
        // Right now it doesn't do shit.
    }

    // called at the end of a report step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::
    endReportStep() {
        // Right now it doesn't do shit.
    }

    // Get the last report step
    template<typename TypeTag>
    const SimulatorReport& 
    BlackoilAquiferModel<TypeTag>::
    lastReport() const {
        return last_report_;
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::
    assemble(const int iterationIdx,
             const double dt)
    {
        // Still doesn't do shit
    }

    // Protected function: Update the primary variables
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::
    updatePrimaryVariables() {
        // Right now it doesn't do shit.
    }

    // Protected function: Init the primary variables
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::
    initPrimaryVariablesEvaluation() const {
        // Right now it doesn't do shit.
    }

    // Protected function: Return number of components in the model.
    template<typename TypeTag>
    int
    BlackoilAquiferModel<TypeTag>::numComponents() const
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
        // Not implemented yet!!!!!!!!!!
        return -1;
        // return aquifers()->number_of_aquifers;
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
    BlackoilAquiferModel<TypeTag>::
    flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
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
    BlackoilAquiferModel<TypeTag>::
    assembleAquiferEq(const double dt)
    {
        // Still doesn't do shit
    }

    // Protected function
    // some preparation work, mostly related to group control and RESV,
    // at the beginning of each time step (Not report step)
    template<typename TypeTag>
    void BlackoilAquiferModel<TypeTag>::
    prepareTimeStep()
    {
        // Not implemented yet!
    }

    // Protected function
    /// return true if aquifers are available in the reservoir
    template<typename TypeTag>
    bool
    BlackoilAquiferModel<TypeTag>::
    aquifersActive() const
    {
        return aquifers_active_;
    }

    // Protected function
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::
    setAquifersActive(const bool aquifers_active)
    {
        aquifers_active_ = aquifers_active;
    }

    // Protected function
    /// return true if Aquifers are available on this process
    template<typename TypeTag>
    bool
    BlackoilAquiferModel<TypeTag>::
    localAquifersActive() const
    {
        return numAquifers() > 0;
    }

} // namespace Opm