import os
import pathlib


test_plot_dir = './Rad1D_test_plots'

data_dir = pathlib.Path(__file__).parents[1]
data_dir = str(data_dir.joinpath('data'))


def pytest_configure(config):
    pass


def pytest_sessionstart(session):
    os.makedirs(test_plot_dir, exist_ok=True)


def pytest_sessionfinish(session, exitstatus):
    pass


def pytest_unconfigure(config):
    pass
