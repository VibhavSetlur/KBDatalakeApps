import os
import argparse
from pathlib import Path
import json
from berdl.query.query_pangenome_local import QueryPangenomeLocal
from berdl.query.query_genome_local import QueryGenomeLocal
from berdl.pangenome.pangenome import BERDLPangenome
from berdl.pangenome.paths_pangenome import PathsPangenome


def main(input_params, selected_clade_member):
    scratch = input_params['_config']['scratch']
    root_pangenome = Path(scratch) / 'pangenome' / selected_clade_member
    paths = PathsPangenome(root=root_pangenome)
    print(paths.root)
    query_pg = QueryPangenomeLocal('/data/reference_data/berdl_db/ke-pangenomes')
    query_g = QueryGenomeLocal()
    berld_pangenome = BERDLPangenome(query_pg, query_g, paths)
    berld_pangenome.run(selected_clade_member)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run BERDL pangenome pipeline"
    )
    parser.add_argument(
        "input_params",
        help="Path to input params JSON file"
    )
    parser.add_argument(
        "selected_clade_member",
        help="clade member to build pangenome"
    )
    #  read input params
    args = parser.parse_args()
    filename_input_params = args.input_params
    selected_clade_member = args.selected_clade_member

    if not os.path.exists(filename_input_params):
        raise FileNotFoundError(
            f"Input params file not found: {filename_input_params}"
        )

    with open(filename_input_params, 'r') as fh:
        input_params = json.load(fh)

    main(input_params, selected_clade_member)
