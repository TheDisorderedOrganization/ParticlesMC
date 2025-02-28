using Arianna
using ParticlesMC

include("parse_input.jl")

function main(args)

    verbose = args["verbose"]
    # Print arguments
    if verbose
        println("Parsed args:")
        for (arg, val) in args
            println("  $arg  =>  $val")
        end
    end

    # Load init file (or files)
    init_path = args["init_file"]
    chains = load_init_files(init_path; args=args, verbose=verbose, tails=".xyz")

    return nothing
    
end

if abspath(PROGRAM_FILE) == @__FILE__
    cl_args = parse_commandline()
    args = load_config_file(cl_args)
    main(args)
end