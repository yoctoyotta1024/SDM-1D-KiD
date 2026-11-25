"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: generate_inputfiles.py
Project: pysdm
Created Date: Monday 23rd November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script to generate a real configuration file (by copying and editing a draft one),
which PySDM requires as model input data.
"""

# %% ---------------------------- GLOBAL IMPORTS ----------------------------- #
import logging

logger = logging.getLogger(__name__)


# %% ------------------------- FUNCTION DEFINITIONS -------------------------- #
def parse_arguments():
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config_filename", type=Path, help="Path to destination (real) YAML config file"
    )
    parser.add_argument(
        "src_config_filename",
        type=Path,
        help="Path to source (draft) YAML config file",
    )
    parser.add_argument(
        "config_params",
        type=json.loads,
        help="Dictionary of configuration parameters to set",
    )
    return parser.parse_args()


def generate_configuration_file(src_config_filename, config_filename, config_params):
    import shutil

    from ruamel.yaml import YAML

    ### --- copy src_config_filename into config_filename --- ##
    config_filename.unlink(missing_ok=True)  # delete any existing config
    shutil.copy(src_config_filename, config_filename)

    ### --- edit selected parameters of config_filename --- ###
    def set_value(node, key, value):
        """Set first 'key' found (at any depth) in YAML 'node' to 'value'."""
        if isinstance(node, dict):
            if key in node:
                node[key] = value
                return True
            return any(set_value(v, key, value) for v in node.values())
        return False

    yaml = YAML()
    with open(config_filename) as file:
        config = yaml.load(file)

    for key, value in config_params.items():
        if not set_value(config, key, value):
            raise ValueError(f"'{key}' not found in {config_filename}")

    with open(config_filename, "w") as file:
        yaml.dump(config, file)


# %% -------------------------------- MAIN ----------------------------------- #
def main(
    config_filename,
    src_config_filename,
    config_params,
):
    if src_config_filename is None or config_filename is None:
        raise ValueError(
            "Cannot generate config file: source and/or destination config filename(s) missing"
        )
    elif not src_config_filename.exists():
        raise ValueError(
            "Cannot generate config file: source config file doesn't exist"
        )
    else:
        logger.info("generating configuration file")
        generate_configuration_file(src_config_filename, config_filename, config_params)


# %% --------------------------- RUN PROGRAM --------------------------------- #
if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.config_filename.resolve(),
        args.src_config_filename.resolve(),
        args.config_params,
    )
