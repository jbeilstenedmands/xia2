from __future__ import annotations

import copy
import logging
import os

import numpy as np

from dials.algorithms.scaling.algorithm import ScalingAlgorithm
from dials.algorithms.scaling.scaling_library import determine_best_unit_cell
from dials.array_family import flex
from dxtbx.model import ExperimentList

logger = logging.getLogger("dials")

from dials.util.options import ArgumentParser
from libtbx import phil


class BatchScale(object):
    def __init__(self, experiments, reflections, reference):
        # super().__init__(events=["run_cosym", "performed_unit_cell_clustering"])
        phil_scope = phil.parse(
            """
    include scope dials.command_line.scale.phil_scope
""",
            process_includes=True,
        )
        parser = ArgumentParser(phil=phil_scope, check_format=False)
        params, _ = parser.parse_args(args=[], quick_parse=True)
        params.model = "KB"
        self.params = params
        self.input_sg = copy.deepcopy(experiments[0][0].crystal.get_space_group())
        self.reference = reference
        self.input_experiments = experiments
        self.input_reflections = reflections
        self._experiments = ExperimentList([])
        self._reflections = []
        self._output_expt_files = []
        self._output_refl_files = []
        self._new_table = copy.deepcopy(self.input_reflections)
        self._new_expts = copy.deepcopy(self.input_experiments)

        # create a single table and expt per batch
        all_expts = ExperimentList([])
        for expts in experiments:
            all_expts.extend(expts)
        best_unit_cell = determine_best_unit_cell(all_expts)
        # self._experiments = all_expts

        for i, (table, expts) in enumerate(zip(self._new_table, self._new_expts)):
            wavelength = np.mean([expt.beam.get_wavelength() for expt in expts])
            expt = copy.deepcopy(expts[0])
            expt.beam.set_wavelength(wavelength)
            expt.crystal.set_unit_cell(best_unit_cell)
            expt.scaling_model = None
            expt.identifier = str(i)
            self._experiments.append(expt)
            # make a new table, apply existing scales

            table["intensity.sum.value.original"] = table["intensity.sum.value"]
            table["intensity.sum.variance.original"] = table["intensity.sum.variance"]
            if "inverse_scale_factor" in table:
                table["inverse_scale_factor_original"] = table["inverse_scale_factor"]
                table["inverse_scale_factor_variance_original"] = table[
                    "inverse_scale_factor_variance"
                ]

                table["intensity.sum.value"] /= table["inverse_scale_factor"]
                table["intensity.sum.variance"] /= (table["inverse_scale_factor"]) ** 2
                del table["inverse_scale_factor"]
                del table["inverse_scale_factor_variance"]
            for k in list(table.experiment_identifiers().keys()):
                del table.experiment_identifiers()[k]
            table["id"] = flex.int(table.size(), i)
            table.experiment_identifiers()[i] = str(i)
            self._reflections.append(table)

        self.algorithm = ScalingAlgorithm(params, self._experiments, self._reflections)

    def run(self):
        self.algorithm.run()
        if self.reference:
            from dials.algorithms.symmetry.reindex_to_reference import (
                determine_reindex_operator_against_reference,
            )
            from dials.util.reference import intensities_from_reference_file

            wavelength = np.mean(
                [expt.beam.get_wavelength() for expt in self._experiments]
            )

            reference_miller_set = intensities_from_reference_file(
                os.fspath(self.reference), wavelength=wavelength
            )
            test_miller_set = self.algorithm.scaled_miller_array
            change_of_basis_op = determine_reindex_operator_against_reference(
                test_miller_set, reference_miller_set
            )
            for i, (expts, refl) in enumerate(
                zip(self.input_experiments, self.input_reflections)
            ):
                for expt in expts:
                    expt.crystal = expt.crystal.change_basis(change_of_basis_op)
                    expt.crystal.set_space_group(self.input_sg)
                expts.as_file(f"processed_{i}.expt")

                refl["miller_index"] = change_of_basis_op.apply(refl["miller_index"])
                refl.as_file(f"processed_{i}.refl")
                self._output_expt_files.append(f"processed_{i}.expt")
                self._output_refl_files.append(f"processed_{i}.refl")
