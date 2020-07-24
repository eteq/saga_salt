import saga_salt
from matplotlib import pyplot as plt

def load_and_plot(fn):
    spec = saga_salt.load_salt_ascii_spec(fn)
    saga_salt.plot_spectrum(spec)
    return spec


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('file')

    args = parser.parse_args()

    load_and_plot(args.file)
    plt.show()