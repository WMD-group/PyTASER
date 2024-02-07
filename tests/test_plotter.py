"""
Tests for the plotting functionality of PyTASer - implicitly tests much of the parsing/generation
functionality as well.
"""

import os

import pytest

from pytaser import plotter
from pytaser.plotter import ev_to_lambda, lambda_to_ev

_file_path = os.path.dirname(__file__)
_DATA_DIR = os.path.join(_file_path, "data_gaas")
_CDTE_DATA_DIR = os.path.join(_file_path, "data_cdte")
_TIO2_DAS_DATA_DIR = os.path.join(_file_path, "data_tio2_das")


def test_ev_to_lambda():
    input_ev = 2.5
    assert round(ev_to_lambda(input_ev), 2) == 495.94


def test_lambda_to_ev():
    input_lambda = 495.94
    assert round(lambda_to_ev(input_lambda), 2) == 2.5


def test_cutoff_transitions(plotter_gaas):
    highest_transitions = [(-2, 1), (-1, 1), (0, 1)]
    highest_transitions.sort()
    relevant_transitions = plotter.cutoff_transitions(
        plotter_gaas.jdos_diff_if, cutoff=0.75, ind_xmin=0, ind_xmax=-1
    )
    relevant_transitions = [x for x in relevant_transitions if x is not None]
    relevant_transitions.sort()
    assert relevant_transitions == highest_transitions


def test_get_plot(plotter_gaas):
    assert plotter_gaas.bandgap_lambda == plotter.ev_to_lambda(plotter_gaas.bandgap)
    assert plotter_gaas.energy_mesh_lambda.all() == plotter.ev_to_lambda(plotter_gaas.energy_mesh_ev).all()


# to run the following image comparison tests and see relative differences, use this in CLI:
# <pytest --mpl --mpl-generate-summary=html test_plotter.py>

# to generate a new set of images for the remote baseline plots, use this in CLI:
# <pytest --mpl-generate-path=data_gaas/remote_baseline_plots test_plotter.py>
# <pytest --mpl-generate-path=data_cdte/remote_baseline_plots test_plotter.py>
# <pytest --mpl-generate-path=data_tio2_das/remote_baseline_plots test_plotter.py>


@pytest.mark.mpl_image_compare(
    style="default",
    baseline_dir=f"{_TIO2_DAS_DATA_DIR}/remote_baseline_plots",
    filename="TiO2_das_energy.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight", "dpi": 100},
)
def test_get_plot_das_ev(plotter_tio2_das):
    """Test get_plot() DAS function for TiO2 with a electronvolts xaxis."""
    return plotter_tio2_das.get_plot(
        xaxis="energy",
        transition_cutoff=0.01,
        xmin=0,
        xmax=10,
        yaxis="das",
    )


@pytest.mark.mpl_image_compare(
    style="default",
    baseline_dir=f"{_TIO2_DAS_DATA_DIR}/remote_baseline_plots",
    filename="TiO2_das_nm.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight", "dpi": 100},
)
def test_get_plot_das_nm(plotter_tio2_das):
    """Test get_plot() DAS function for TiO2 with a wavelength xaxis."""
    return plotter_tio2_das.get_plot(
        xaxis="wavelength",
        transition_cutoff=0.01,
        xmin=200,
        xmax=1200,
        yaxis="das",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_DATA_DIR}/remote_baseline_plots",
    filename="tas_ev_gaas.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_tas_ev(plotter_gaas):
    """
    Test get_plot() TAS function for GaAs with a 25% cutoff and electronvolts xaxis.
    """
    return plotter_gaas.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        transition_cutoff=0.75,
        xmin=None,
        xmax=None,
        ymin=None,
        ymax=None,
        yaxis="tas",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_DATA_DIR}/remote_baseline_plots",
    filename="tas_lambda_gaas.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_tas_lambda(plotter_gaas):
    """Test get_plot() TAS function for GaAs with a 25% cutoff and a wavelength xaxis."""
    return plotter_gaas.get_plot(
        relevant_transitions="auto",
        xaxis="wavelength",
        transition_cutoff=0.75,
        xmin=None,
        xmax=1200,
        ymin=None,
        ymax=None,
        yaxis="tas",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_DATA_DIR}/remote_baseline_plots",
    filename="jdos_ev_gaas.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_jdos_ev(plotter_gaas):
    """Test get_plot() JDOS function for GaAs with a 25% cutoff and a electronvolts xaxis."""
    return plotter_gaas.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        transition_cutoff=0.75,
        xmin=None,
        xmax=None,
        ymin=None,
        ymax=None,
        yaxis="jdos",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_DATA_DIR}/remote_baseline_plots",
    filename="jdos_lambda_gaas.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_jdos_lambda(plotter_gaas):
    """Test get_plot() JDOS function for GaAs with a 25% cutoff and a wavelength xaxis."""
    return plotter_gaas.get_plot(
        relevant_transitions="auto",
        xaxis="wavelength",
        transition_cutoff=0.75,
        xmin=None,
        xmax=1200,
        ymin=None,
        ymax=None,
        yaxis="jdos",
    )


## Test with CdTe, with many more transitions:
def test_get_plot_cdte(plotter_cdte):
    assert plotter_cdte.bandgap_lambda == plotter.ev_to_lambda(plotter_cdte.bandgap)
    assert plotter_cdte.energy_mesh_lambda.all() == plotter.ev_to_lambda(plotter_cdte.energy_mesh_ev).all()


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="tas_ev_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_tas_ev_cdte(plotter_cdte):
    """Test get_plot() TAS function for CdTe with the default cutoff and electronvolts xaxis."""
    return plotter_cdte.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        xmin=0,
        xmax=5,
        yaxis="tas",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="tas_lambda_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_tas_lambda_cdte(plotter_cdte):
    """Test get_plot() TAS function for CdTe with the default cutoff and wavelength xaxis."""
    return plotter_cdte.get_plot(
        relevant_transitions="auto",
        xaxis="wavelength",
        xmin=None,
        xmax=1200,
        yaxis="tas",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="jdos_ev_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_jdos_ev_cdte(plotter_cdte):
    """Test get_plot() JDOS function for CdTe with the default cutoff and electronvolts xaxis."""
    return plotter_cdte.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        xmin=0,
        xmax=5,
        yaxis="jdos",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="jdos_lambda_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_jdos_lambda_cdte(plotter_cdte):
    """Test get_plot() JDOS function for CdTe with the default cutoff and wavelength xaxis."""
    return plotter_cdte.get_plot(
        relevant_transitions="auto",
        xaxis="wavelength",
        xmin=None,
        xmax=1200,
        yaxis="jdos",
    )


def test_line_color_consistency(plotter_cdte):
    """Test that the same transition has the same color in all plots, when transition_cutoff is
    changed.
    """
    fig = plotter_cdte.get_plot()  # with default transition_cutoff of 0.03
    line = [_line for _line in fig.gca().lines if "(-2, 0)" in _line.get_label()][0]
    line_color = line.get_color()

    fig = plotter_cdte.get_plot(transition_cutoff=0.3)  # this removes lines before (-2, 0)
    line = [_line for _line in fig.gca().lines if "(-2, 0)" in _line.get_label()][0]
    assert line_color == line.get_color()

    # check for JDOS plots as well:
    fig = plotter_cdte.get_plot(yaxis="jdos")  # with default transition_cutoff of 0.03
    line = [_line for _line in fig.gca().lines if "(-2, 0)" in _line.get_label()][0]
    line_color = line.get_color()

    fig = plotter_cdte.get_plot(transition_cutoff=0.3, yaxis="jdos")  # removes lines before (-2, 0)
    line = [_line for _line in fig.gca().lines if "(-2, 0)" in _line.get_label()][0]
    assert line_color == line.get_color()


def test_get_plot_alpha_cdte_no_waveder(plotter_cdte):
    """Test informative error for get_plot() yaxis="alpha" when no WAVEDER was parsed."""
    with pytest.raises(
        ValueError,
        match="The `alpha` option for yaxis can only be chosen if the "
        "TASGenerator object was created using VASP outputs!",
    ):
        plotter_cdte.get_plot(
            relevant_transitions="auto",
            xaxis="energy",
            xmin=0,
            xmax=5,
            yaxis="alpha",
        )


def test_get_plot_tas_absorption_only_cdte_no_waveder(plotter_cdte):
    """Test informative error for get_plot() yaxis="tas_absorption_only" when no WAVEDER was parsed."""
    with pytest.raises(
        ValueError,
        match="The `tas_absorption_only` option for yaxis can only be chosen if the "
        "TASGenerator object was created using VASP outputs!",
    ):
        plotter_cdte.get_plot(
            relevant_transitions="auto",
            xaxis="energy",
            xmin=0,
            xmax=5,
            yaxis="tas_absorption_only",
        )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="tas_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_tas_cdte(plotter_cdte_vasp):
    """Test get_plot() yaxis="tas" for CdTe with the default cutoff and electronvolts xaxis."""
    return plotter_cdte_vasp.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        xmin=0,
        xmax=5,
        yaxis="tas",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="jdos_diff_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_jdos_diff_cdte_vasp(plotter_cdte_vasp):
    """
    Test get_plot() yaxis="jdos_diff" function for CdTe with the default cutoff and electronvolts xaxis.
    """
    return plotter_cdte_vasp.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        xmin=0,
        xmax=5,
        yaxis="jdos_diff",
    )


# jdos_diff plot should be the same for both plotter_cdte_vasp and plotter_cdte_vasp_vr_only, so compare
# to same baseline:
@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="jdos_diff_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_jdos_diff_cdte_vasp_vr_only(plotter_cdte_vasp_vr_only):
    """Test get_plot() yaxis="jdos_diff" for CdTe with the default cutoff and electronvolts xaxis."""
    return plotter_cdte_vasp_vr_only.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        xmin=0,
        xmax=5,
        yaxis="jdos_diff",
    )


# jdos_diff plot should be the same as "tas" for plotter_cdte_vasp_vr_only (because no WAVEDER parsed),
# so compare to same baseline:
@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="jdos_diff_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_jdos_diff_cdte(plotter_cdte_vasp_vr_only):
    """Test get_plot() yaxis="tas" for CdTe with the default cutoff and electronvolts xaxis."""
    return plotter_cdte_vasp_vr_only.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        xmin=0,
        xmax=5,
        yaxis="tas",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="alpha_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_alpha_cdte(plotter_cdte_vasp):
    """Test get_plot() yaxis="alpha" for CdTe with the default cutoff and electronvolts xaxis."""
    return plotter_cdte_vasp.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        xmin=0,
        xmax=5,
        yaxis="alpha",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="tas_absorption_only_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_tas_absorption_only_cdte_vasp(plotter_cdte_vasp):
    """
    Test get_plot() yaxis="tas_absorption_only" function for CdTe
    with the default cutoff and electronvolts xaxis.
    """
    return plotter_cdte_vasp.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        xmin=0,
        xmax=5,
        yaxis="tas_absorption_only",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="jdos_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_jdos_cdte(plotter_cdte_vasp):
    """Test get_plot() yaxis="jdos" function for CdTe with the default cutoff and electronvolts xaxis."""
    return plotter_cdte_vasp.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        xmin=0,
        xmax=5,
        yaxis="jdos",
    )


# jdos should be the same for plotter_cdte_vasp and plotter_cdte_vasp_vr_only, so compare to same baseline:
@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="jdos_cdte.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_jdos_cdte_vasp_vr_only(plotter_cdte_vasp_vr_only):
    """Test get_plot() yaxis="jdos" for CdTe with the default cutoff and electronvolts xaxis."""
    return plotter_cdte_vasp_vr_only.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        xmin=0,
        xmax=5,
        yaxis="jdos",
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="tas_cdte_custom_legend.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_tas_cdte_custom_legend(plotter_cdte_vasp):
    """Test get_plot() yaxis="tas" for CdTe with kwargs for plt.legend()."""
    return plotter_cdte_vasp.get_plot(
        relevant_transitions="auto",
        xaxis="energy",
        xmin=0,
        xmax=5,
        yaxis="tas",
        # kwargs for plt.legend():
        ncols=2,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.0),
        frameon=False,
        borderaxespad=0.0,
        handlelength=1.5,
        handletextpad=0.5,
    )


@pytest.mark.mpl_image_compare(
    baseline_dir=f"{_CDTE_DATA_DIR}/remote_baseline_plots",
    filename="tas_cdte_scissor_bandgap_3eV.png",
    savefig_kwargs={"transparent": True, "bbox_inches": "tight"},
)
def test_get_plot_tas_cdte_bandgap_scissor(plotter_cdte_bg_3_vasp):
    """Test get_plot() for CdTe with scissored bandgap to 3 eV."""
    return plotter_cdte_bg_3_vasp.get_plot(xaxis="energy")
