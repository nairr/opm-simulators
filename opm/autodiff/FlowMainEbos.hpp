/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
  Copyright 2014 STATOIL ASA.

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

#ifndef OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED
#define OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED

#include <opm/simulators/ParallelFileMerger.hpp>
#include <opm/simulators/ensureDirectoryExists.hpp>

#include <opm/autodiff/BlackoilModelEbos.hpp>
#include <opm/autodiff/NewtonIterationBlackoilSimple.hpp>
#include <opm/autodiff/NewtonIterationBlackoilCPR.hpp>
#include <opm/autodiff/NewtonIterationBlackoilInterleaved.hpp>
#include <opm/autodiff/MissingFeatures.hpp>
#include <opm/autodiff/moduleVersion.hpp>
#include <opm/autodiff/ExtractParallelGridInformationToISTL.hpp>
#include <opm/autodiff/RedistributeDataHandles.hpp>

#include <opm/core/props/satfunc/RelpermDiagnostics.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>
#include <opm/common/ResetLocale.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>

#include <ewoms/version.hh>

namespace Opm
{

    /// \brief Gather cell data to global random access iterator
    /// \tparam ConstIter The type of constant iterator.
    /// \tparam Iter The type of the mutable iterator.
    /// \param grid The distributed CpGrid where loadbalance has been run.
    /// \param local The local container from which the data should be sent.
    /// \param global The global container to gather to.
    /// \warning The global container has to have the correct size!
    template<class ConstIter, class Iter>
    void gatherCellDataToGlobalIterator(const Dune::CpGrid& grid,
                                         const ConstIter& local_begin,
                                         const Iter& global_begin)
    {
#if HAVE_MPI
        FixedSizeIterCopyHandle<ConstIter,Iter> handle(local_begin,
                                                   global_begin);
        const auto& gatherScatterInf = grid.cellScatterGatherInterface();
        Dune::VariableSizeCommunicator<> comm(grid.comm(),
                                              gatherScatterInf);
        comm.backward(handle);
#endif
    }


    // The FlowMain class is the ebos based black-oil simulator.
    class FlowMainEbos
    {
    public:
        typedef TTAG(EclFlowProblem) TypeTag;
        typedef typename GET_PROP(TypeTag, MaterialLaw)::EclMaterialLawManager MaterialLawManager;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator) EbosSimulator;
        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

        typedef Opm::SimulatorFullyImplicitBlackoilEbos Simulator;
        typedef typename Simulator::ReservoirState ReservoirState;
        typedef typename Simulator::OutputWriter OutputWriter;

        /// This is the main function of Flow.
        /// It runs a complete simulation, with the given grid and
        /// simulator classes, based on user command-line input.  The
        /// content of this function used to be in the main() function of
        /// flow.cpp.
        int execute(int argc, char** argv)
        {
            try {
                // we always want to use the default locale, and thus spare us the trouble
                // with incorrect locale settings.
                resetLocale();

                setupParallelism(argc, argv);
                printStartupMessage();
                const bool ok = setupParameters(argc, argv);
                if (!ok) {
                    return EXIT_FAILURE;
                }

                setupEbosSimulator();
                setupOutput();
                setupLogging();
                extractMessages();
                setupGridAndProps();
                runDiagnostics();
                setupState();
                writeInit();
                setupOutputWriter();
                setupLinearSolver();
                createSimulator();

                // Run.
                auto ret =  runSimulator();

                mergeParallelLogFiles();

                return ret;
            }
            catch (const std::exception &e) {
                std::ostringstream message;
                message  << "Program threw an exception: " << e.what();

                if( output_cout_ )
                {
                    // in some cases exceptions are thrown before the logging system is set
                    // up.
                    if (OpmLog::hasBackend("STREAMLOG")) {
                        OpmLog::error(message.str());
                    }
                    else {
                        std::cout << message.str() << "\n";
                    }
                }

                return EXIT_FAILURE;
            }
        }

    protected:
        void setupParallelism(int argc, char** argv)
        {
            // MPI setup.
            // Must ensure an instance of the helper is created to initialise MPI.
            // For a build without MPI the Dune::FakeMPIHelper is used, so rank will
            // be 0 and size 1.
            const Dune::MPIHelper& mpi_helper = Dune::MPIHelper::instance(argc, argv);
            mpi_rank_ = mpi_helper.rank();
            const int mpi_size = mpi_helper.size();
            output_cout_ = ( mpi_rank_ == 0 );
            must_distribute_ = ( mpi_size > 1 );

#ifdef _OPENMP
            // OpenMP setup.
            if (!getenv("OMP_NUM_THREADS")) {
                // Default to at most 4 threads, regardless of
                // number of cores (unless ENV(OMP_NUM_THREADS) is defined)
                int num_cores = omp_get_num_procs();
                int num_threads = std::min(4, num_cores);
                omp_set_num_threads(num_threads);
            }
#pragma omp parallel
            if (omp_get_thread_num() == 0) {
                // omp_get_num_threads() only works as expected within a parallel region.
                const int num_omp_threads = omp_get_num_threads();
                if (mpi_size == 1) {
                    std::cout << "OpenMP using " << num_omp_threads << " threads." << std::endl;
                } else {
                    std::cout << "OpenMP using " << num_omp_threads << " threads on MPI rank " << mpi_rank_ << "." << std::endl;
                }
            }
#endif
        }

        // Print startup message if on output rank.
        void printStartupMessage()
        {
            if (output_cout_) {
                const int lineLen = 70;
                const std::string version = moduleVersionName();
                const std::string banner = "This is flow (version "+version+")";
                const int bannerPreLen = (lineLen - 2 - banner.size())/2;
                const int bannerPostLen = bannerPreLen + (lineLen - 2 - banner.size())%2;
                std::cout << "**********************************************************************\n";
                std::cout << "*                                                                    *\n";
                std::cout << "*" << std::string(bannerPreLen, ' ') << banner << std::string(bannerPostLen, ' ') << "*\n";
                std::cout << "*                                                                    *\n";
                std::cout << "* Flow is a simulator for fully implicit three-phase black-oil flow, *\n";
                std::cout << "*            and is part of OPM. For more information see:           *\n";
                std::cout << "*                       http://opm-project.org                       *\n";
                std::cout << "*                                                                    *\n";
                std::cout << "**********************************************************************\n\n";
            }
        }

        // Read parameters, see if a deck was specified on the command line, and if
        // it was, insert it into parameters.
        // Writes to:
        //   param_
        // Returns true if ok, false if not.
        bool setupParameters(int argc, char** argv)
        {
            param_ = ParameterGroup(argc, argv, false, output_cout_);

            // See if a deck was specified on the command line.
            if (!param_.unhandledArguments().empty()) {
                if (param_.unhandledArguments().size() != 1) {
                    std::cerr << "You can only specify a single input deck on the command line.\n";
                    return false;
                } else {
                    const auto casename = this->simulationCaseName( param_.unhandledArguments()[ 0 ] );
                    param_.insertParameter("deck_filename", casename.string() );
                }
            }

            // We must have an input deck. Grid and props will be read from that.
            if (!param_.has("deck_filename")) {
                std::cerr << "This program must be run with an input deck.\n"
                    "Specify the deck filename either\n"
                    "    a) as a command line argument by itself\n"
                    "    b) as a command line parameter with the syntax deck_filename=<path to your deck>, or\n"
                    "    c) as a parameter in a parameter file (.param or .xml) passed to the program.\n";
                return false;
            }
            return true;
        }

        // Set output_to_files_ and set/create output dir. Write parameter file.
        // Writes to:
        //   output_to_files_
        //   output_dir_
        // Throws std::runtime_error if failed to create (if requested) output dir.
        void setupOutput()
        {
            output_to_files_ = output_cout_ && param_.getDefault("output", true);

            // Setup output directory.
            auto& ioConfig = eclState().getIOConfig();
            // Default output directory is the directory where the deck is found.
            const std::string default_output_dir = ioConfig.getOutputDir();
            output_dir_ = param_.getDefault("output_dir", default_output_dir);
            // Override output directory if user specified.
            ioConfig.setOutputDir(output_dir_);

            // Write parameters used for later reference. (only if rank is zero)
            if (output_to_files_) {
                // Create output directory if needed.
                ensureDirectoryExists(output_dir_);
                // Write simulation parameters.
                param_.writeParam(output_dir_ + "/simulation.param");
            }
        }

        // Setup OpmLog backend with output_dir.
        void setupLogging()
        {
            std::string deck_filename = param_.get<std::string>("deck_filename");
            // create logFile
            using boost::filesystem::path;
            path fpath(deck_filename);
            std::string baseName;
            std::ostringstream debugFileStream;
            std::ostringstream logFileStream;

            if (boost::to_upper_copy(path(fpath.extension()).string()) == ".DATA") {
                baseName = path(fpath.stem()).string();
            } else {
                baseName = path(fpath.filename()).string();
            }

            logFileStream << output_dir_ << "/" << baseName;
            debugFileStream << output_dir_ << "/" << "." << baseName;

            if ( must_distribute_ && mpi_rank_ != 0 )
            {
                // Added rank to log file for non-zero ranks.
                // This prevents message loss.
                debugFileStream << "."<< mpi_rank_;
                // If the following file appears then there is a bug.
                logFileStream << "." << mpi_rank_;
            }
            logFileStream << ".PRT";
            debugFileStream << ".DEBUG";

            std::string debugFile = debugFileStream.str();
            logFile_ = logFileStream.str();

            std::shared_ptr<EclipsePRTLog> prtLog = std::make_shared<EclipsePRTLog>(logFile_ , Log::NoDebugMessageTypes, false, output_cout_);
            std::shared_ptr<StreamLog> streamLog = std::make_shared<StreamLog>(std::cout, Log::StdoutMessageTypes);
            OpmLog::addBackend( "ECLIPSEPRTLOG" , prtLog );
            OpmLog::addBackend( "STREAMLOG", streamLog);
            std::shared_ptr<StreamLog> debugLog = std::make_shared<EclipsePRTLog>(debugFile, Log::DefaultMessageTypes, false, output_cout_);
            OpmLog::addBackend( "DEBUGLOG" ,  debugLog);
            const auto& msgLimits = eclState().getSchedule().getMessageLimits();
            const std::map<int64_t, int> limits = {{Log::MessageType::Note, msgLimits.getCommentPrintLimit(0)},
                                                   {Log::MessageType::Info, msgLimits.getMessagePrintLimit(0)},
                                                   {Log::MessageType::Warning, msgLimits.getWarningPrintLimit(0)},
                                                   {Log::MessageType::Error, msgLimits.getErrorPrintLimit(0)},
                                                   {Log::MessageType::Problem, msgLimits.getProblemPrintLimit(0)},
                                                   {Log::MessageType::Bug, msgLimits.getBugPrintLimit(0)},
                                                   {Log::MessageType::Probleminfo, msgLimits.getProbleminfoPrintLimit(0)},
                                                   {Log::MessageType::Warninginfo, msgLimits.getWarninginfoPrintLimit(0)}};
            prtLog->setMessageLimiter(std::make_shared<MessageLimiter>());
            prtLog->setMessageFormatter(std::make_shared<SimpleMessageFormatter>(false));
            streamLog->setMessageLimiter(std::make_shared<MessageLimiter>(10, limits));
            streamLog->setMessageFormatter(std::make_shared<SimpleMessageFormatter>(true));

            // Read parameters.
            if ( output_cout_ )
            {
                OpmLog::debug("\n---------------    Reading parameters     ---------------\n");
            }
        }

        void mergeParallelLogFiles()
        {
            // force closing of all log files.
            OpmLog::removeAllBackends();

            if( mpi_rank_ != 0 || !must_distribute_ || !output_to_files_ )
            {
                return;
            }

            namespace fs = boost::filesystem;
            fs::path output_path(".");
            if ( param_.has("output_dir") )
            {
                output_path = fs::path(output_dir_);
            }

            fs::path deck_filename(param_.get<std::string>("deck_filename"));

            std::for_each(fs::directory_iterator(output_path),
                          fs::directory_iterator(),
                          detail::ParallelFileMerger(output_path, deck_filename.stem().string()));
        }

        void setupEbosSimulator()
        {
            std::string progName("flow_ebos");
            std::string deckFile("--ecl-deck-file-name=");
            deckFile += param_.get<std::string>("deck_filename");
            char* ptr[2];
            ptr[ 0 ] = const_cast< char * > (progName.c_str());
            ptr[ 1 ] = const_cast< char * > (deckFile.c_str());
            EbosSimulator::registerParameters();
            Ewoms::setupParameters_< TypeTag > ( 2, ptr );
            ebosSimulator_.reset(new EbosSimulator(/*verbose=*/false));
            ebosSimulator_->model().applyInitialSolution();

            // Create a grid with a global view.
            globalGrid_.reset(new Grid(grid()));
            globalGrid_->switchToGlobalView();

            try {
                if (output_cout_) {
                    MissingFeatures::checkKeywords(deck());
                }

                // Possible to force initialization only behavior (NOSIM).
                if (param_.has("nosim")) {
                    const bool nosim = param_.get<bool>("nosim");
                    auto& ioConfig = eclState().getIOConfig();
                    ioConfig.overrideNOSIM( nosim );
                }
            }
            catch (const std::invalid_argument& e) {
                std::cerr << "Failed to create valid EclipseState object. See logfile: " << logFile_ << std::endl;
                std::cerr << "Exception caught: " << e.what() << std::endl;
                throw;
            }

            // Possibly override IOConfig setting (from deck) for how often RESTART files should get written to disk (every N report step)
            if (param_.has("output_interval")) {
                const int output_interval = param_.get<int>("output_interval");
                eclState().getRestartConfig().overrideRestartWriteInterval( size_t( output_interval ) );
            }
        }

        // Create distributed property objects.
        // Writes to:
        //   fluidprops_
        void setupGridAndProps()
        {
            Dune::CpGrid& grid = ebosSimulator_->gridManager().grid();

            // create the legacy properties objects
            fluidprops_.reset(new BlackoilPropsAdFromDeck(deck(),
                                                          eclState(),
                                                          materialLawManager(),
                                                          grid));

            // Geological properties
            bool use_local_perm = param_.getDefault("use_local_perm", true);
            geoprops_.reset(new DerivedGeology(grid, *fluidprops_, eclState(), use_local_perm, &ebosProblem().gravity()[0]));
        }

        const Deck& deck() const
        { return ebosSimulator_->gridManager().deck(); }

        Deck& deck()
        { return ebosSimulator_->gridManager().deck(); }

        const EclipseState& eclState() const
        { return ebosSimulator_->gridManager().eclState(); }

        EclipseState& eclState()
        { return ebosSimulator_->gridManager().eclState(); }

        // Initialise the reservoir state. Updated fluid props for SWATINIT.
        // Writes to:
        //   state_
        //   threshold_pressures_
        //   fluidprops_ (if SWATINIT is used)
        void setupState()
        {
            const PhaseUsage pu = Opm::phaseUsageFromDeck(deck());
            const Grid& grid = this->grid();

            // Need old-style fluid object for init purposes (only).
            BlackoilPropertiesFromDeck props(deck(),
                                             eclState(),
                                             materialLawManager(),
                                             grid.size(/*codim=*/0),
                                             grid.globalCell().data(),
                                             grid.logicalCartesianSize().data(),
                                             param_);


            // Init state variables (saturation and pressure).
            if (param_.has("init_saturation")) {
                state_.reset(new ReservoirState(grid.size(/*codim=*/0),
                                                grid.numFaces(),
                                                props.numPhases()));

                initStateBasic(grid.size(/*codim=*/0),
                               grid.globalCell().data(),
                               grid.logicalCartesianSize().data(),
                               grid.numFaces(),
                               Opm::UgGridHelpers::faceCells(grid),
                               Opm::UgGridHelpers::beginFaceCentroids(grid),
                               Opm::UgGridHelpers::beginCellCentroids(grid),
                               Grid::dimension,
                               props, param_, gravity(), *state_);

                initBlackoilSurfvol(Opm::UgGridHelpers::numCells(grid), props, *state_);

                enum { Oil = BlackoilPhases::Liquid, Gas = BlackoilPhases::Vapour };
                if (pu.phase_used[Oil] && pu.phase_used[Gas]) {
                    const int numPhases = props.numPhases();
                    const int numCells  = Opm::UgGridHelpers::numCells(grid);

                    // Uglyness 1: The state is a templated type, here we however make explicit use BlackoilState.
                    auto& gor = state_->getCellData( BlackoilState::GASOILRATIO );
                    const auto& surface_vol = state_->getCellData( BlackoilState::SURFACEVOL );
                    for (int c = 0; c < numCells; ++c) {
                        // Uglyness 2: Here we explicitly use the layout of the saturation in the surface_vol field.
                        gor[c] = surface_vol[ c * numPhases + pu.phase_pos[Gas]] / surface_vol[ c * numPhases + pu.phase_pos[Oil]];
                    }
                }
            } else if (deck().hasKeyword("EQUIL")) {
                // Which state class are we really using - what a f... mess?
                state_.reset( new ReservoirState( Opm::UgGridHelpers::numCells(grid),
                                                  Opm::UgGridHelpers::numFaces(grid),
                                                  props.numPhases()));

                initStateEquil(grid, props, deck(), eclState(), gravity(), *state_);
                //state_.faceflux().resize(Opm::UgGridHelpers::numFaces(grid), 0.0);
            } else {
                state_.reset( new ReservoirState( Opm::UgGridHelpers::numCells(grid),
                                                  Opm::UgGridHelpers::numFaces(grid),
                                                  props.numPhases()));
                initBlackoilStateFromDeck(Opm::UgGridHelpers::numCells(grid),
                                          Opm::UgGridHelpers::globalCell(grid),
                                          Opm::UgGridHelpers::numFaces(grid),
                                          Opm::UgGridHelpers::faceCells(grid),
                                          Opm::UgGridHelpers::beginFaceCentroids(grid),
                                          Opm::UgGridHelpers::beginCellCentroids(grid),
                                          Opm::UgGridHelpers::dimensions(grid),
                                          props, deck(), gravity(), *state_);
            }

            // The capillary pressure is scaled in fluidprops_ to match the scaled capillary pressure in props.
            if (deck().hasKeyword("SWATINIT")) {
                const int numCells = Opm::UgGridHelpers::numCells(grid);
                std::vector<int> cells(numCells);
                for (int c = 0; c < numCells; ++c) { cells[c] = c; }
                std::vector<double> pc = state_->saturation();
                props.capPress(numCells, state_->saturation().data(), cells.data(), pc.data(), nullptr);
                fluidprops_->setSwatInitScaling(state_->saturation(), pc);
            }
            initHydroCarbonState(*state_, pu, Opm::UgGridHelpers::numCells(grid), deck().hasKeyword("DISGAS"), deck().hasKeyword("VAPOIL"));
        }

        // Extract messages from parser.
        // Writes to:
        //    OpmLog singleton.
        void extractMessages()
        {
            if ( !output_cout_ )
            {
                return;
            }

            auto extractMessage = [this](const Message& msg) {
                auto log_type = this->convertMessageType(msg.mtype);
                const auto& location = msg.location;
                if (location) {
                    OpmLog::addMessage(log_type, Log::fileMessage(location.filename, location.lineno, msg.message));
                } else {
                    OpmLog::addMessage(log_type, msg.message);
                }
            };

            // Extract messages from Deck.
            for(const auto& msg : deck().getMessageContainer()) {
                extractMessage(msg);
            }

            // Extract messages from EclipseState.
            for (const auto& msg : eclState().getMessageContainer()) {
                extractMessage(msg);
            }
        }

        // Run diagnostics.
        // Writes to:
        //   OpmLog singleton.
        void runDiagnostics()
        {
            if( ! output_cout_ )
            {
                return;
            }

            // Run relperm diagnostics
            RelpermDiagnostics diagnostic;
            diagnostic.diagnosis(eclState(), deck(), this->grid());
        }

        void writeInit()
        {
            bool output      = param_.getDefault("output", true);
            bool output_ecl  = param_.getDefault("output_ecl", true);
            if( output && output_ecl )
            {
                const Grid& grid = this->globalGrid();

                if( output_cout_ ){
                    const EclipseGrid& inputGrid = eclState().getInputGrid();
                    eclIO_.reset(new EclipseIO(eclState(), UgGridHelpers::createEclipseGrid( grid , inputGrid )));
                }

                const NNC* nnc = &geoprops_->nonCartesianConnections();
                data::Solution globaltrans;

                if ( must_distribute_ )
                {
                    // dirty and dangerous hack!
                    // We rely on opmfil in GeoProps being hardcoded to true
                    // which prevents the pinch processing from running.
                    // Ergo the nncs are unchanged.
                    nnc = &eclState().getInputNNC();

                    // Gather the global simProps
                    data::Solution localtrans = geoprops_->simProps(this->grid());
                    for( const auto& localkeyval: localtrans)
                    {
                        auto& globalval = globaltrans[localkeyval.first].data;
                        const auto& localval  = localkeyval.second.data;

                        if( output_cout_ )
                        {
                            globalval.resize( grid.size(0));
                        }
                        gatherCellDataToGlobalIterator(this->grid(), localval.begin(),
                                                       globalval.begin());
                    }
                }
                else
                {
                    globaltrans = geoprops_->simProps(grid);
                }

                if( output_cout_ )
                {
                eclIO_->writeInitial(globaltrans,
                                              *nnc);
                }
            }
        }

        // Setup output writer.
        // Writes to:
        //   output_writer_
        void setupOutputWriter()
        {
            // create output writer after grid is distributed, otherwise the parallel output
            // won't work correctly since we need to create a mapping from the distributed to
            // the global view
            output_writer_.reset(new OutputWriter(grid(),
                                                  param_,
                                                  eclState(),
                                                  std::move(eclIO_),
                                                  Opm::phaseUsageFromDeck(deck())) );
        }

        // Run the simulator.
        // Returns EXIT_SUCCESS if it does not throw.
        int runSimulator()
        {
            const auto& schedule = eclState().getSchedule();
            const auto& timeMap = schedule.getTimeMap();
            auto& ioConfig = eclState().getIOConfig();
            SimulatorTimer simtimer;

            // initialize variables
            const auto& initConfig = eclState().getInitConfig();
            simtimer.init(timeMap, (size_t)initConfig.getRestartStep());

            if (!ioConfig.initOnly()) {
                if (output_cout_) {
                    std::string msg;
                    msg = "\n\n================ Starting main simulation loop ===============\n";
                    OpmLog::info(msg);
                }

                SimulatorReport successReport = simulator_->run(simtimer, *state_);
                SimulatorReport failureReport = simulator_->failureReport();

                if (output_cout_) {
                    std::ostringstream ss;
                    ss << "\n\n================    End of simulation     ===============\n\n";
                    successReport.reportFullyImplicit(ss, &failureReport);
                    OpmLog::info(ss.str());
                    if (param_.anyUnused()) {
                        // This allows a user to catch typos and misunderstandings in the
                        // use of simulator parameters.
                        std::cout << "--------------------   Unused parameters:   --------------------\n";
                        param_.displayUsage();
                        std::cout << "----------------------------------------------------------------" << std::endl;
                    }
                }

                if (output_to_files_) {
                    std::string filename = output_dir_ + "/walltime.txt";
                    std::fstream tot_os(filename.c_str(), std::fstream::trunc | std::fstream::out);
                    successReport.reportParam(tot_os);
                }
            } else {
                if (output_cout_) {
                    std::cout << "\n\n================ Simulation turned off ===============\n" << std::flush;
                }

            }
            return EXIT_SUCCESS;
        }

        // Setup linear solver.
        // Writes to:
        //   fis_solver_
        void setupLinearSolver()
        {
            typedef typename BlackoilModelEbos :: ISTLSolverType ISTLSolverType;

            extractParallelGridInformationToISTL(grid(), parallel_information_);
            fis_solver_.reset( new ISTLSolverType( param_, parallel_information_ ) );
        }

        /// This is the main function of Flow.
        // Create simulator instance.
        // Writes to:
        //   simulator_
        void createSimulator()
        {
            // Create the simulator instance.
            simulator_.reset(new Simulator(*ebosSimulator_,
                                           param_,
                                           *geoprops_,
                                           *fluidprops_,
                                           *fis_solver_,
                                           FluidSystem::enableDissolvedGas(),
                                           FluidSystem::enableVaporizedOil(),
                                           eclState(),
                                           *output_writer_,
                                           defunctWellNames()));
        }

    private:
        boost::filesystem::path simulationCaseName( const std::string& casename ) {
            namespace fs = boost::filesystem;

            const auto exists = []( const fs::path& f ) -> bool {
                if( !fs::exists( f ) ) return false;

                if( fs::is_regular_file( f ) ) return true;

                return fs::is_symlink( f )
                && fs::is_regular_file( fs::read_symlink( f ) );
            };

            auto simcase = fs::path( casename );

            if( exists( simcase ) ) {
                return simcase;
            }

            for( const auto& ext : { std::string("data"), std::string("DATA") } ) {
                if( exists( simcase.replace_extension( ext ) ) ) {
                    return simcase;
                }
            }

            throw std::invalid_argument( "Cannot find input case " + casename );
        }


        int64_t convertMessageType(const Message::type& mtype)
        {
            switch (mtype) {
            case Message::type::Debug:
                return Log::MessageType::Debug;
            case Message::type::Info:
                return Log::MessageType::Info;
            case Message::type::Warning:
                return Log::MessageType::Warning;
            case Message::type::Error:
                return Log::MessageType::Error;
            case Message::type::Problem:
                return Log::MessageType::Problem;
            case Message::type::Bug:
                return Log::MessageType::Bug;
            case Message::type::Note:
                return Log::MessageType::Note;
            case Message::type::Probleminfo:
                return Log::MessageType::Probleminfo;
            case Message::type::Warninginfo:
                return Log::MessageType::Warninginfo;
            }
            throw std::logic_error("Invalid messages type!\n");
        }

        Grid& grid()
        { return ebosSimulator_->gridManager().grid(); }

        const Grid& globalGrid()
        { return *globalGrid_; }

        Problem& ebosProblem()
        { return ebosSimulator_->problem(); }

        const Problem& ebosProblem() const
        { return ebosSimulator_->problem(); }

        std::shared_ptr<MaterialLawManager> materialLawManager()
        { return ebosProblem().materialLawManager(); }

        Scalar gravity() const
        { return ebosProblem().gravity()[2]; }

        std::unordered_set<std::string> defunctWellNames() const
        { return ebosSimulator_->gridManager().defunctWellNames(); }

        std::unique_ptr<EbosSimulator> ebosSimulator_;
        int  mpi_rank_ = 0;
        bool output_cout_ = false;
        bool must_distribute_ = false;
        ParameterGroup param_;
        bool output_to_files_ = false;
        std::string output_dir_ = std::string(".");
        std::unique_ptr<BlackoilPropsAdFromDeck> fluidprops_;
        std::unique_ptr<DerivedGeology> geoprops_;
        std::unique_ptr<ReservoirState> state_;
        std::unique_ptr<EclipseIO> eclIO_;
        std::unique_ptr<OutputWriter> output_writer_;
        boost::any parallel_information_;
        std::unique_ptr<NewtonIterationBlackoilInterface> fis_solver_;
        std::unique_ptr<Simulator> simulator_;
        std::string logFile_;
        // Needs to be shared pointer because it gets initialzed before MPI_Init.
        std::shared_ptr<Grid> globalGrid_;
    };
} // namespace Opm

#endif // OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED
