import argparse
import runpy
import sys


class CLICommand:
    """Run GPAW's Python interpreter."""

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('--dry-run', type=int, default=0,
                            metavar='NCPUS',
                            help='Dry run on NCPUS cpus.')
        group = parser.add_mutually_exclusive_group()
        group.add_argument('--command', '-c', dest='cmd',
                           help='Program passed in as string.')
        group.add_argument('--module', '-m',
                           help='Run library module as a script.')
        parser.add_argument('arguments', metavar='ARG',
                            help='Arguments passed to program in '
                            'sys.argv[1:].',
                            nargs=argparse.REMAINDER)

    @staticmethod
    def run(args):
        if args.cmd:
            sys.argv[:] = ['-c'] + args.arguments
            d = {}
            exec(args.cmd, d, d)
        elif args.module:
            sys.argv[:] = [args.module] + args.arguments
            runpy.run_module(args.module, run_name='__main__')
        else:
            sys.argv[:] = args.arguments
            runpy.run_path(args.arguments[0], run_name='__main__')
