class CLICommand:
    short_description = "Run GPAW's Python interpreter"

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('arguments', nargs='*')

    @staticmethod
    def run(args):
        pass