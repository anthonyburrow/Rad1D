import os


test_plot_dir = './Rad1D_test_plots'


def pytest_configure(config):
    pass


def pytest_sessionstart(session):
    os.makedirs(test_plot_dir, exist_ok=True)


def pytest_sessionfinish(session, exitstatus):
    pass


def pytest_unconfigure(config):
    pass
