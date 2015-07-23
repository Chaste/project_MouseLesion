/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTFIBROBLASTS1D_HPP_
#define TESTFIBROBLASTS1D_HPP_

/*
 * = Two dimensional mouse ventricle with lesion simulation =
 *
 * This is the code that was used to perform the simulation in Mahoney ''et al.'' (2015).
 *
 * == Code Walkthrough ==
 */

#include <cxxtest/TestSuite.h>

/*
 *  These two includes are in this project, and not a standard part of Chaste v3.3
 */
#include "ScarCellFactory.hpp"
#include "ScarConductivityModifier.hpp"

/*
 * The rest of these includes are standard Chaste files...
 */
#include "GmshMeshReader.hpp"
#include "MonodomainProblem.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "CellProperties.hpp" // For analysing APs

#include "PetscSetupAndFinalize.hpp"

class TestFibroblasts1d : public CxxTest::TestSuite
{
private:
    /** Length of the 1D domain we're simulating. */
    double mRegionWidth;

    /** Length of the scar in the centre. */
    double mScarRadius;

    /** Width of the boundary region (myocytes but different conductivity) */
    double mBoundaryWidth;

public:
    void Test1dPropagation() throw (Exception)
    {
#ifdef CHASTE_CVODE
        if (*(CommandLineArguments::Instance()->p_argc) == 1)
        {
            std::cout << "TestFibroblasts takes the following arguments:\n"
                    " Mandatory:\n"
                    " * --scar-conduction-scaling <x> The proportion of the intra-cellular conduction to\n"
                    "                                 apply in the central scar region.\n"
                    "\n"
                    " Optional:\n"
                    " * --use-neutral-cells <true or false>  Whether to use neutral cells OR,\n"
                    " * --use-fibroblasts <true or false>    a fibroblast electrophysiology model,\n"
                    "                                        in the scar region.\n"
                    " * --scar-membrane-area-scaling <x> the scaling factor for membrane area (per unit vol) in scar region\n"
                    " * --boundary-conduction-scaling <x> The proportion of intra-cellular conduction to\n"
                    "                                     apply in the boundary of the scar (still myocytes).\n"
                    "\n"
                    " * --pacing-period <x>  The time between paces applied on x=0 (defaults to 100ms)\n"
                    " * --end-time <x>       How long to perform the simulation for (defaults to pacing period).\n"
                    " * --lesion-pacing      Whether to perform pacing in the centre of the lesion (default false).\n";
            return;
        }

        /**
         * SET OPTIONS FOR THIS RUN
         *
         * These are the defaults.
         */
        double pacing_cycle_length = 100;
        bool use_neutral_cell_model = false;
        bool use_fibroblasts = false;
        bool lesion_pacing = false; // Alternative is remote/sinus pacing.
        ScarShape shape = CIRCLE;
        double scar_membrane_area_scaling = 1.0;
        const std::string output_folder = "FibroblastSims1d/";

        if (CommandLineArguments::Instance()->OptionExists("--pacing-period"))
        {
            pacing_cycle_length = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--pacing-period");
        }
        double end_time = pacing_cycle_length; // Default end time is pacing cycle length

        if (CommandLineArguments::Instance()->OptionExists("--end-time"))
        {
            end_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--end-time");
        }

        // Here we say what proportion of the normal conductivity we want in the scar
        double scaling_factor;
        if (CommandLineArguments::Instance()->OptionExists("--scar-conduction-scaling"))
        {
            scaling_factor = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--scar-conduction-scaling");
            assert(scaling_factor >= 0.0);
            assert(scaling_factor <= 1.0);
        }
        else
        {
            EXCEPTION("Please provide the command line option '--scar-conduction-scaling' with a value between 0 and 1.");
        }

        double scaling_factor_boundary = 1.0;
        if (CommandLineArguments::Instance()->OptionExists("--boundary-conduction-scaling"))
        {
            scaling_factor_boundary = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--boundary-conduction-scaling");
            assert(scaling_factor_boundary >= 0.0);
            assert(scaling_factor_boundary <= 1.0);
        }

        // And decide whether to use a neutral cell model, or normal mytocytes.
        if (CommandLineArguments::Instance()->OptionExists("--use-neutral-cells"))
        {
            use_neutral_cell_model = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("--use-neutral-cells");
        }
        // And decide whether to use a fibroblast model
        if (CommandLineArguments::Instance()->OptionExists("--use-fibroblasts"))
        {
            use_fibroblasts = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("--use-fibroblasts");
        }
        if (use_neutral_cell_model && use_fibroblasts)
        {
            EXCEPTION("You can only set '--use-neutral-cells true' OR '--use-fibroblasts true', not both.");
        }

        if (CommandLineArguments::Instance()->OptionExists("--scar-membrane-area-scaling"))
        {
            scar_membrane_area_scaling = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--scar-membrane-area-scaling");
            assert(scar_membrane_area_scaling >= 0.0);
            assert(scar_membrane_area_scaling <= 1.0);
        }

        if (CommandLineArguments::Instance()->OptionExists("--lesion-pacing"))
        {
            lesion_pacing = true;
        }

        // Set up a unique output folder for these options
        std::stringstream sub_directory;
        if (use_neutral_cell_model)
        {
            sub_directory << "Neutral";
        }
        else if (use_fibroblasts)
        {
            sub_directory << "Fibroblast";
        }
        else
        {
            sub_directory << "Myocyte";
        }

        sub_directory << "_scar_cond_" << scaling_factor << "_cap_" << scar_membrane_area_scaling << "_boundary_cond_" << scaling_factor_boundary << "_period_" << pacing_cycle_length << "_end_" << end_time;

        if (lesion_pacing)
        {
            sub_directory << "_LesionPacing";
        }

        /*
         * SET UP MESH
         */
        DistributedTetrahedralMesh<1,1> mesh;

        // Mesh properties we will use in the cell factory and conductivity modifier.
        mRegionWidth = 0.5; //cm
        mScarRadius = 0.1; //cm
        mBoundaryWidth = 0.01; //cm

        double h=0.005;
        mesh.ConstructRegularSlabMesh(h, mRegionWidth, mRegionWidth);

        /*
         * SET STANDARD OPTIONS
         */
        HeartConfig::Instance()->SetSimulationDuration(end_time); //ms
        HeartConfig::Instance()->SetOutputDirectory(output_folder + sub_directory.str());
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);

        /*
         * NUMERICAL METHOD PARAMETERS
         */
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 0.1, 0.1);

        /*
         * SET UP HETEROGENEOUS CONDUCTIVITY
         */
        double intra_conductivity = 1.75; // Chaste Defaults
        double extra_conductivity = 7.0; // Chaste Defaults

        // Set up the defaults.
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(intra_conductivity, intra_conductivity));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(extra_conductivity, extra_conductivity));

        /*
         * COMPLETE THE SET UP OF PROBLEM
         */

        // Make a scar cell factory
        ScarCellFactory<1> cell_factory(mRegionWidth,
                                        mScarRadius,
                                        pacing_cycle_length,
                                        scar_membrane_area_scaling,  // chi scaling factor.
                                        use_neutral_cell_model,
                                        use_fibroblasts,
                                        shape, // doesn't make any difference in 1D
                                        lesion_pacing);

        // Create a conductivity modifier
        ScarConductivityModifier<1> modifier(&mesh,
                                             shape, // doesn't make any difference in 1D
                                             mRegionWidth,
                                             mScarRadius,
                                             scaling_factor, // D scaling factor.
                                             scar_membrane_area_scaling,  // chi scaling factor.
                                             mBoundaryWidth,
                                             scaling_factor_boundary, // D scaling factor (chi not implemented for boundary yet).
                                             true, // whether a lesion applied
                                             false);

        MonodomainProblem<1> problem( &cell_factory );
        problem.SetMesh( &mesh );
        
        problem.SetWriteInfo(); // Writes max and min voltages out.
        problem.Initialise();
        problem.GetTissue()->SetConductivityModifier(&modifier);
        problem.Solve();

        // Print out a couple of node traces to file for easy plotting.
        // Output file handler should be called collectively, the 'false' says "don't wipe the folder!"
        OutputFileHandler handler(output_folder + sub_directory.str(), false);

        // These methods must be called in parallel.
        std::vector<unsigned> node_indices;
        node_indices.push_back(cell_factory.GetNodeIndexCentreOfScar()); // 0
        node_indices.push_back(cell_factory.GetNodeIndexLeftTissue()); // 2
        node_indices.push_back(cell_factory.GetNodeIndexRightTissue()); // 3
        node_indices.push_back(cell_factory.GetNodeIndexLeftOfScar()); // 4
        node_indices.push_back(cell_factory.GetNodeIndexRightOfScar()); // 5

        // But we only want the master to write out the actual file.
        if (PetscTools::AmMaster())
        {
            FileFinder hdf5_dir(output_folder + sub_directory.str(), RelativeTo::ChasteTestOutput);
            Hdf5DataReader reader(hdf5_dir, "results");

            // Retrieve the voltage traces at nodes of interest
            std::vector<std::vector<double> > voltage_traces;
            for (unsigned i=0; i<node_indices.size(); i++)
            {
                voltage_traces.push_back(reader.GetVariableOverTime("V", node_indices[i]));
            }

            // Output the full voltage traces at these nodes
            std::vector<double> times = reader.GetUnlimitedDimensionValues();
            out_stream output_file = handler.OpenOutputFile("voltage_traces_at_selected_nodes.dat");
            
            for (unsigned i=0; i<times.size(); i++)
            {
                *output_file << times[i];
                for (unsigned j=0; j<node_indices.size(); j++)
                {
                    *output_file << "\t" << voltage_traces[j][i] ;
                }
                *output_file << "\n";
            }
            output_file->close();

            // Also output some of our action potential summary statistics
            output_file = handler.OpenOutputFile("action_potential_properties_at_selected_nodes.dat");
            *output_file << "Node\tAPD90(ms)\tAPD70(ms)\tAPD50(ms)\tMaxUpstroke(mV/ms)\tPeakVoltage(mV)\tAmplitude(mV)\tResting(mV)" << std::endl;
            const double voltage_AP_threshold = -70.0;
            for (unsigned i=0; i<node_indices.size(); i++)
            {
                CellProperties voltage_properties(voltage_traces[i], times, voltage_AP_threshold);
                double apd90, apd70, apd50, upstroke, peak, amplitude, resting;
                try
                {
                    apd90 = voltage_properties.GetLastActionPotentialDuration(90);
                    apd70 = voltage_properties.GetLastActionPotentialDuration(70);
                    apd50 = voltage_properties.GetLastActionPotentialDuration(50);
                    upstroke = voltage_properties.GetLastCompleteMaxUpstrokeVelocity();
                    peak = voltage_properties.GetLastCompletePeakPotential();
                    amplitude = voltage_properties.GetLastActionPotentialAmplitude();
                    resting = voltage_properties.GetRestingPotentials().back();
                }
                catch (Exception& e)
                {   // In case there is no action potential at all.
                    apd90 = 0.0;
                    apd70 = 0.0;
                    apd50 = 0.0;
                    upstroke = 0.0;
                    peak = 0.0;
                    amplitude = 0.0;
                    resting = voltage_traces[i].back();
                }
                *output_file << i << "\t"
                             << apd90 << "\t" << apd70 << "\t" << apd50 << "\t"
                             << upstroke << "\t" << peak << "\t" << amplitude << "\t" << resting << std::endl;
            }
            output_file->close();
        }

#else
        std::cout << "CVODE is not installed, or CHASTE is not configured to use it, check your hostconfig settings." << std::endl;
        // TS_ASSERT(false); // uncomment if you want to ensure CVODE is set up on your system.
#endif // CHASTE_CVODE
    }
};

#endif // TESTFIBROBLASTS_HPP_
