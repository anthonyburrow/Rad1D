import os


test_plot_dir = './Rad1D_tests'


def pytest_sessionstart(session):
    os.makedirs(test_plot_dir, exist_ok=True)
