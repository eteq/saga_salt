import saga_salt

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('file')

    args = parser.parse_args()

    spec = saga_salt.load_salt_ascii_spec(args.file)
    saga_salt.plot_spectrum(spec)