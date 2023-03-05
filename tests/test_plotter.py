import pytest
from pytaser import plotter
from pytaser.plotter import ev_to_lambda, lambda_to_ev
import os

_file_path = os.path.dirname(__file__)
_DATA_DIR = os.path.join(_file_path, "data_gaas")


def test_ev_to_lambda():
    input_ev = 2.5
    assert round(ev_to_lambda(input_ev), 2) == 495.94


def test_lamda_to_ev():
    input_lamda = 495.94
    assert round(lambda_to_ev(input_lamda), 2) == 2.5


def test_cutoff_transitions(plotter_gaas):
    highest_transitions = [(-2, 1), (-1, 1), (0, 1)]
    highest_transitions.sort()
    relevant_transitions = plotter.cutoff_transitions(plotter_gaas.tas_decomp, cutoff=0.75, ind_xmin=0, ind_xmax=-1)
    relevant_transitions.sort()
    assert relevant_transitions == highest_transitions


def test_get_plot(plotter_gaas):
    assert plotter_gaas.bandgap_lambda == plotter.ev_to_lambda(plotter_gaas.bandgap_ev)
    assert plotter_gaas.energy_mesh_lambda.all() == plotter.ev_to_lambda(plotter_gaas.energy_mesh_ev).all()

#to run the following image comparison tests and see relative differences, use the CLI
# "<pytest --mpl --mpl-generate-summary=html test_plotter.py>"

@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_DATA_DIR}/remote_baseline_plots",
    filename="tas_ev_gaas.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_tas_ev(plotter_gaas):
    """Test get_plot() TAS function for GaAs with a 25% cutoff and a electronvolts xaxis"""
    fig = plotter_gaas.get_plot(relevant_transitions="auto",
                                xaxis="energy",
                                transition_cutoff=0.75,
                                xmin=None,
                                xmax=None,
                                ymin=None,
                                ymax=None,
                                yaxis="tas"
                                )
    return fig


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_DATA_DIR}/remote_baseline_plots",
    filename="tas_lambda_gaas.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_tas_lambda(plotter_gaas):
    """Test get_plot() TAS function for GaAs with a 25% cutoff and a wavelength xaxis"""
    fig = plotter_gaas.get_plot(relevant_transitions="auto",
                                xaxis="wavelength",
                                transition_cutoff=0.75,
                                xmin=None,
                                xmax=1200,
                                ymin=None,
                                ymax=None,
                                yaxis="tas"
                                )
    return fig


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_DATA_DIR}/remote_baseline_plots",
    filename="jdos_ev_gaas.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_jdos_ev(plotter_gaas):
    """Test get_plot() JDOS function for GaAs with a 25% cutoff and a electronvolts xaxis"""
    fig = plotter_gaas.get_plot(relevant_transitions="auto",
                                xaxis="energy",
                                transition_cutoff=0.75,
                                xmin=None,
                                xmax=None,
                                ymin=None,
                                ymax=None,
                                yaxis="jdos"
                                )
    return fig


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_DATA_DIR}/remote_baseline_plots",
    filename="jdos_lambda_gaas.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_jdos_lambda(plotter_gaas):
    """Test get_plot() JDOS function for GaAs with a 25% cutoff and a wavelength xaxis"""
    fig = plotter_gaas.get_plot(relevant_transitions="auto",
                                xaxis="wavelength",
                                transition_cutoff=0.75,
                                xmin=None,
                                xmax=1200,
                                ymin=None,
                                ymax=None,
                                yaxis="jdos"
                                )
    return fig
