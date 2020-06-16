from . import arguments
from . import assembly_stats
from . import draft_model
from . import util


def entry():
    # Get command line arguments
    args = arguments.parse()

    # Execute workflows
    if args.command is None:
        pass
    elif args.command == 'assembly_qc':
        print(assembly_stats.run(args.assembly_fp))
    elif args.command == 'annotation':
        pass
    elif args.command == 'draft_model':
        model = util.read_model_and_check(args.ref_model_fp, args.ref_gbk_fp)
        draft_model.run(args.assembly_fp, args.ref_gbk_fp, model, args.output_fp)
    elif args.command == 'model_fba':
        pass
    else:
        assert False


if __name__ == '__main__':
    entry()
