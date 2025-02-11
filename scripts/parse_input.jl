using ArgParse
using YAML

function parse_commandline()
    parser = ArgParseSettings()
    @add_arg_table! parser begin
        "init_file"
        help = "Path to the initial configuration file (accepts multiple files)"
        arg_type = String
        required = true
        "config_file"
        help = "Path to the simulation parameters file (e.g., input.yaml)"
        arg_type = String
        required = true
        "--steps"
        help = "Ovveride the number of steps in the config file"
        arg_type = Int
        "--nsim"
        help = "Ovveride the number of chains per config file"
        arg_type = Int
        "--temperature", "-T"
        help = "Ovveride the temperature in the input file"
        arg_type = Float64
        "--density", "-D"
        help = "Ovveride the density in the input file (affine transformation)"
        arg_type = Float64
        "--model"
        help = "Ovveride the model in the input file"
        arg_type = String
        "--list_type"
        help = "Ovveride the cell list type (EmptyList or LinkedList)"
        arg_type = String
        "--verbose", "-v"
        help = "verbose"
        action = :store_true
        "--seed"
        help = "Override random number seed"
        arg_type = Int
    end

    return parse_args(parser)

end

function load_config_file(cl_args)
    args = YAML.load_file(cl_args["config_file"])
    for key in keys(cl_args)
        if cl_args[key] !== nothing
            args[key] = cl_args[key]
        end
    end
    return args
end

nothing