import os
from pathlib import Path
import json
from berdl.genome_paths import GenomePaths
from berdl.prep_genome_set import BERDLPreGenome
from cobrakbase import KBaseAPI


def main():
    #  read input params from scratch
    # input_params =

    #  setup clients/methods
    kbase = KBaseAPI('')

    #  extract genomes
    genomes = []

    paths = GenomePaths(root=Path("/test/syncom").resolve())
    berdl_prep_genomes = BERDLPreGenome(kbase, paths)
    user_to_clade, ani_clade, df_ani_fitness, df_ani_phenotype = berdl_prep_genomes.run(genomes)


if __name__ == "__main__":
    main()
