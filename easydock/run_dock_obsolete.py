import sys
import easydock.run_dock as run_dock


def main():
    sys.stderr.write('Warning: "run_dock" is an obsolete command, please use "easydock" instead. "run_dock" will be '
                     'removed in future versions.\n\n')
    run_dock.main()


if __name__ == '__main__':
    main()
