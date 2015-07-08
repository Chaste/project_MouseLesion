#!/bin/bash

cd ../..
CHASTE_PATH=`pwd`

# Run the VTK time annotations script on all these
# and copy the results back to where we can find them and play with them.
if test -d ${CHASTE_TEST_OUTPUT}/FibroblastSims; then
    cd ${CHASTE_TEST_OUTPUT}/FibroblastSims
    # Copy the results to a sensible place for the analysis script
    RESULTS_DIR=${CHASTE_PATH}/projects/MouseLesion/test/data/results_capacitance
    mkdir -p ${RESULTS_DIR}
    for dir in *;
      do (cd $dir && echo $dir && ${CHASTE_PATH}/python/utils/AddVtuTimeAnnotations.py vtk_output/results.vtu vtk_output/annotated.vtu && cp voltage_traces_at_selected_nodes.dat ${RESULTS_DIR}/$dir.dat && cp action_potential_properties_at_selected_nodes.dat ${RESULTS_DIR}/${dir}_APs.dat && cd ..);
    done
fi

if test -d ${CHASTE_TEST_OUTPUT}/FibroblastSims3d; then
    cd ${CHASTE_TEST_OUTPUT}/FibroblastSims3d
    for dir in *;
      do (cd $dir && echo $dir && ${CHASTE_PATH}/python/utils/AddVtuTimeAnnotations.py vtk_output/results.vtu vtk_output/annotated.vtu &&  cd ..);
    done
fi

