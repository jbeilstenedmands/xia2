from __future__ import annotations

import contextlib
import copy
import json
import logging
import os
from pathlib import Path
from typing import List, Tuple

import iotbx.phil
from cctbx import crystal, sgtbx, uctbx
from dials.algorithms.indexing.ssx.analysis import (
    generate_html_report,
    generate_plots,
    make_summary_table,
    report_on_crystal_clusters,
)
from dials.algorithms.integration.ssx.ssx_integrate import (
    generate_html_report as generate_integration_html_report,
)
from dials.algorithms.shoebox import MaskCode
from dials.array_family import flex
from dials.command_line.find_spots import working_phil as find_spots_phil
from dials.command_line.ssx_index import index, phil_scope
from dials.command_line.ssx_integrate import run_integration, working_phil
from dials.util import log, tabulate
from dials.util.ascii_art import spot_counts_per_image_plot
from dxtbx.model import ExperimentList
from dxtbx.serialize import load
from xfel.clustering.cluster import Cluster

from xia2.Handlers.Streams import banner

logger = logging.getLogger(__name__)


@contextlib.contextmanager
def run_in_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)


@contextlib.contextmanager
def log_to_file(filename):
    try:
        log.config(logfile=filename, add_console=False)
        dials_logger = logging.getLogger("dials")
        yield dials_logger
    finally:
        dials_logger = logging.getLogger("dials")
        for h in dials_logger.handlers:
            dials_logger.removeHandler(h)


def ssx_find_spots(working_directory: Path) -> flex.reflection_table:

    with run_in_directory(working_directory):
        # Set up the input
        imported_expts = load.experiment_list("imported.expt", check_format=True)
        params = find_spots_phil.extract()
        # Do spot-finding
        logger.notice(banner("Spotfinding"))
        reflections = flex.reflection_table.from_observations(imported_expts, params)
        good = MaskCode.Foreground | MaskCode.Valid
        reflections["n_signal"] = reflections["shoebox"].count_mask_values(good)
        logger.info(spot_counts_per_image_plot(reflections))
    return reflections


def ssx_index(
    working_directory: Path,
    nproc: int = 1,
    space_group: sgtbx.space_group = None,
    unit_cell: uctbx.unit_cell = None,
) -> Tuple[ExperimentList, flex.reflection_table, List[Cluster]]:

    with run_in_directory(working_directory):
        logfile = "dials.ssx_index.log"
        with log_to_file(logfile) as dials_logger:
            # Set up the input and log it to the dials log file
            strong_refl = flex.reflection_table.from_file("strong.refl")
            imported_expts = load.experiment_list("imported.expt", check_format=False)
            params = phil_scope.extract()
            params.indexing.nproc = nproc
            input_ = (
                "Input parameters:\n  reflections = strong.refl\n"
                + "  experiments = imported.expt\n  indexing.nproc = {nproc}\n"
            )
            if unit_cell:
                params.indexing.known_symmetry.unit_cell = unit_cell
                uc = ",".join(str(i) for i in unit_cell.parameters())
                input_ += f"  indexing.known_symmetry.unit_cell = {uc}\n"
            if space_group:
                params.indexing.known_symmetry.space_group = space_group
                input_ += (
                    f"  indexing.known_symmetry.space_group = {str(space_group)}\n"
                )
            dials_logger.info(input_)

            # Do the indexing
            logger.notice(banner("Indexing"))
            indexed_experiments, indexed_reflections, summary_data = index(
                imported_expts, strong_refl, params
            )
            n_images = len({e.imageset.get_path(0) for e in indexed_experiments})
            report = (
                f"{indexed_reflections.size()} spots indexed on {n_images} images\n"
                + "\nSummary of images sucessfully indexed\n"
                + make_summary_table(summary_data)
            )
            logger.info(report)
            dials_logger.info(report)

            # Report on clustering, and generate html report and json output
            crystal_symmetries = [
                crystal.symmetry(
                    unit_cell=expt.crystal.get_unit_cell(),
                    space_group=expt.crystal.get_space_group(),
                )
                for expt in indexed_experiments
            ]
            cluster_plots, large_clusters = report_on_crystal_clusters(
                crystal_symmetries, True
            )
            summary_plots = generate_plots(summary_data)
            summary_plots.update(cluster_plots)
            generate_html_report(summary_plots, "dials.ssx_index.html")
            with open("dials.ssx_index.json", "w") as outfile:
                json.dump(summary_plots, outfile)
    return indexed_experiments, indexed_reflections, large_clusters


def generate_refinement_step_table(refiner):
    import math

    rmsd_multipliers = []
    header = ["Step", "Nref"]
    for (name, units) in zip(refiner._target.rmsd_names, refiner._target.rmsd_units):
        if units == "mm":
            header.append(name + "\n(mm)")
            rmsd_multipliers.append(1.0)
        elif units == "rad":  # convert radians to degrees for reporting
            header.append(name + "\n(deg)")
            rmsd_multipliers.append(180 / math.pi)
        else:  # leave unknown units alone
            header.append(name + "\n(" + units + ")")

    rows = []
    for i in range(refiner._refinery.history.get_nrows()):
        rmsds = [
            r * m
            for (r, m) in zip(refiner._refinery.history["rmsd"][i], rmsd_multipliers)
        ]
        rows.append(
            [str(i), str(refiner._refinery.history["num_reflections"][i])]
            + [f"{r:.5g}" for r in rmsds]
        )
    return tabulate(rows, header)


def run_refinement(working_directory: Path) -> None:
    logger.notice(banner("Joint refinement"))

    from dials.command_line.refine import run_dials_refine, working_phil

    with run_in_directory(working_directory):
        logfile = "dials.refine.log"
        with log_to_file(logfile) as dials_logger:
            params = working_phil.extract()
            indexed_refl = flex.reflection_table.from_file("indexed.refl")
            indexed_expts = load.experiment_list("indexed.expt", check_format=False)
            params.refinement.parameterisation.beam.fix = "all"
            params.refinement.parameterisation.auto_reduction.action = "fix"
            params.refinement.parameterisation.detector.fix_list = "Tau1"
            params.refinement.refinery.engine = "SparseLevMar"
            params.refinement.reflections.outlier.algorithm = "sauter_poon"
            expts, refls, refiner, _ = run_dials_refine(
                indexed_expts, indexed_refl, params
            )
            dials_logger.info("Saving refined experiments to refined.expt")
            expts.as_file("refined.expt")
            dials_logger.info(
                "Saving reflections with updated predictions to refined.refl"
            )
            refls.as_file("refined.refl")
            step_table = generate_refinement_step_table(refiner)
            logger.info("Summary of joint refinement steps:\n" + step_table)


def ssx_integrate(
    working_directory: Path, integration_params: iotbx.phil.scope_extract
) -> List[Cluster]:

    with run_in_directory(working_directory):
        logfile = "dials.ssx_integrate.log"
        with log_to_file(logfile) as dials_logger:
            # Set up the input and log it to the dials log file
            indexed_refl = flex.reflection_table.from_file(
                "indexed.refl"
            ).split_by_experiment_id()
            indexed_expts = load.experiment_list("indexed.expt", check_format=True)
            params = working_phil.extract()
            params.output.batch_size = 100
            params.algorithm = integration_params.algorithm
            input_ = (
                "Input parameters:\n  reflections = indexed.refl\n"
                + f"  experiments = indexed.expt\n  algorithm = {integration_params.algorithm}\n"
                + "  output.batch_size = 100\n"
            )
            if integration_params.algorithm == "ellipsoid":
                model = integration_params.ellipsoid.rlp_mosaicity
                params.profile.ellipsoid.refinement.outlier_probability = 0.95
                params.profile.ellipsoid.refinement.max_separation = 1
                params.profile.ellipsoid.prediction.probability = 0.95
                params.profile.ellipsoid.rlp_mosaicity.model = model
                input_ += "  profile.ellipsoid.refinement.outlier_probability = 0.95\n"
                input_ += "  profile.ellipsoid.refinement.max_separation = 1\n"
                input_ += "  profile.ellipsoid.prediction.probability = 0.95\n"
                input_ += f"  profile.ellipsoid.rlp_mosaicity.model = {model}\n"
            dials_logger.info(input_)

            # Run the integration
            logger.notice(banner("Integrating"))
            integrated_crystal_symmetries = []
            n_refl, n_cryst = (0, 0)
            for i, (int_expt, int_refl, aggregator) in enumerate(
                run_integration(indexed_refl, indexed_expts, params)
            ):
                reflections_filename = f"integrated_{i+1}.refl"
                experiments_filename = f"integrated_{i+1}.expt"
                n_refl += int_refl.size()
                dials_logger.info(
                    f"Saving {int_refl.size()} reflections to {reflections_filename}"
                )
                int_refl.as_file(reflections_filename)
                n_cryst += len(int_expt)
                dials_logger.info(f"Saving the experiments to {experiments_filename}")
                int_expt.as_file(experiments_filename)

                integrated_crystal_symmetries.extend(
                    [
                        crystal.symmetry(
                            unit_cell=copy.deepcopy(cryst.get_unit_cell()),
                            space_group=copy.deepcopy(cryst.get_space_group()),
                        )
                        for cryst in int_expt.crystals()
                    ]
                )
            logger.info(f"{n_refl} reflections integrated from {n_cryst} crystals")

            # Report on clustering, and generate html report and json output
            plots = {}
            cluster_plots, large_clusters = report_on_crystal_clusters(
                integrated_crystal_symmetries,
                make_plots=True,
            )
            plots = aggregator.make_plots()
            plots.update(cluster_plots)
            generate_integration_html_report(plots, "dials.ssx_integrate.html")
            with open("dials.ssx_integrate.json", "w") as outfile:
                json.dump(plots, outfile, indent=2)
    return large_clusters


def best_cell_from_cluster(cluster: Cluster) -> Tuple:
    from cctbx import crystal
    from cctbx.sgtbx.lattice_symmetry import metric_subgroups
    from cctbx.uctbx import unit_cell

    input_symmetry = crystal.symmetry(
        unit_cell=unit_cell(cluster.medians[0:6]), space_group_symbol="P 1"
    )
    group = metric_subgroups(input_symmetry, 3.00).result_groups[0]
    uc_params_conv = group["best_subsym"].unit_cell().parameters()
    sg = group["best_subsym"].space_group_info().symbol_and_number()
    return sg, uc_params_conv


def condensed_unit_cell_info(clusters: List[Cluster]) -> str:
    al, be, ga = "med_" + "\u03B1", "med_" + "\u03B2", "med_" + "\u03B3"
    out_str = f"{'n_xtals'} {'s.g.':>7} {'':>5} {'med_a':>7} {'med_b':>7} {'med_c':>7} {al:>6} {be:>6} {ga:>6}\n"

    for cluster in clusters:
        sorted_pg_comp = sorted(
            cluster.pg_composition.items(),
            key=lambda x: -1 * x[1],
        )
        pg_str = ",".join([str(pg[0]) for pg in sorted_pg_comp])
        p = [f"{i:.2f}" for i in cluster.medians]
        sds = [f"{i:.2f}" for i in cluster.stdevs]
        out_str += f"{len(cluster.members):>7} {pg_str:>7} value {p[0]:>7} {p[1]:>7} {p[2]:>7} {p[3]:>6} {p[4]:>6} {p[5]:>6}\n"
        out_str += f"{'':>7} {'':>7} stdev {sds[0]:>7} {sds[1]:>7} {sds[2]:>7} {sds[3]:>6} {sds[4]:>6} {sds[5]:>6}\n"
    return out_str
