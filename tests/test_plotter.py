import pytest
import numpy as np
from pytaser.kpoints import get_kpoint_weights
from pytaser import plotter
from pytaser.plotter import TASPlotter
from pytaser.tas import Tas
from monty.serialization import loadfn

tas_obj = loadfn('data_sil/example_tas.yaml')
plotter_obj = loadfn('data_sil/example_plotter.yaml')
max_lambda = loadfn('data_sil/max_abs_vals.yaml')
max_val_lambda = loadfn('data_sil/max_val_lambda.json')
cutoffs = loadfn('data_sil/cutoff_list.json')
transition_energies = loadfn('data_sil/relevant_transition_energies.json')


def test_ev_to_lambda():
    input_ev = 2.5
    assert round(ev_to_lambda(input_ev), 2) == 495.94


def test_lamda_to_ev():
    input_lamda = 495.94
    assert round(lambda_to_ev(input_lamda), 2) == 2.5


@pytest.mark.parametrize("tas_object", tas_obj)
@pytest.mark.parametrize("example_plotter_obj", plotter_obj)
@pytest.mark.parametrize("example_max_lambda", max_lambda)
@pytest.mark.parametrize("example_max_val_lambda", max_val_lambda)
@pytest.mark.parametrize("example_cutoffs", cutoffs)
@pytest.mark.parametrize("example_transition_energies", transition_energies)
def test_get_plot(tas_object, example_plotter_obj, example_max_lambda, example_max_val_lambda, example_cutoffs,
                  example_transition_energies):

    # testing the attributes of the TASPlotter object
    plotter_obj = TASPlotter(tas_object)
    assert plotter_obj.bandgap_ev == example_plotter_obj.bandgap_ev
    assert plotter_obj.bandgap_lambda == plotter.ev_to_lambda(plotter_obj.bandgap_ev)
    assert plotter_obj.energy_mesh_lambda == plotter.ev_to_lambda(plotter_obj.energy_mesh_ev)

    # testing the relevant_transitions = "auto" feature in get_plot()
    xmin_lambda = 750
    xmax_lambda = 4000
    xmax_ind_lambda = np.abs(plotter_obj.energy_mesh_lambda - xmin_lambda).argmin()
    xmin_ind_lambda = np.abs(plotter_obj.energy_mesh_lambda - xmax_lambda).argmin()
    max_abs_vals_lambda = {key: np.max(abs(val[xmin_ind_lambda:xmax_ind_lambda])) for key, val in
                           self.tas_decomp.items()}
    max_val_lambda = max(max_abs_vals_lambda.values())
    transition_cutoff = 0.75
    cutoff_list = []
    for transition, value in max_abs_vals_lambda.items():
        if value >= (max_val_lambda * transition_cutoff):
            cutoff_list += [value]
    assert max_abs_vals_lambda == example_max_lambda
    assert max_val_lambda == example_max_val_lambda
    assert len(cutoff_list) == len(example_cutoffs)

    # testing the relevant transitions maps correctly
    transition_array = [(-3, -1), (-1, 0), (0, 1), (1, 3)]
    relevant_transition_energies = []
    for transition in transition_array:
        relevant_transition_energies += [tas_object.tas_decomp[transition][xmin_ind_lambda:xmax_ind_lambda]]
    assert relevant_transition_energies == example_transition_energies
