import os
import sys
import json
import argparse
import traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path = ["/kb/module/lib","/deps/KBUtilLib/src"] + sys.path

from modelseedpy.core.msgenome import MSGenome, MSFeature
import pandas as pd
from KBDatalakeApps.KBDatalakeUtils import KBDataLakeUtils, run_model_reconstruction, run_phenotype_simulation, run_model_reconstruction2


def read_rast_as_genome(filename_rast, genome_id, ontology_term_rast='RAST') -> MSGenome:
    d_rast = pd.read_csv(filename_rast, sep='\t', index_col=0).to_dict()[ontology_term_rast]
    genome = MSGenome()
    genome.id = genome_id
    genome.scientific_name = genome_id
    features = []
    for feature_id, rast_str in d_rast.items():
        feature = MSFeature(feature_id, "")
        if not pd.isna(rast_str):
            for v in rast_str.split('; '):
                feature.add_ontology_term('RAST', v)
        features.append(feature)
    genome.add_features(features)

    return genome


def main(params):
    """Run model reconstruction, phenotype simulation, table building, and BERDL database pipeline."""
    token = params['token']
    scratch = params['scratch']
    kbversion = params['kbversion']
    max_phenotypes = params.get('max_phenotypes', None)
    input_refs = params['input_refs']
    module_path = params.get('module_path', '/kb/module')
    data_path = f"{module_path}/data"
    print('model_pipeline params:', params)

    # Construct SDK config (container paths)
    sdk_config = {
        "kbversion": kbversion,
        "scratch": scratch,
        "reference_path": "/data",
        "module_path": module_path,
        "max_phenotypes": max_phenotypes,
        "experimental_phenotype_datafile": f"{data_path}/experimental_data.json",
        "phenotypeset_file": f"{data_path}/full_phenotype_set.json",
        "fitness_genomes_dir": "/data/reference_data/phenotype_data/fitness_genomes/",
        "reference_phenosim_dir": "/data/reference_data/phenotype_data/phenosims",
        "experimental_essentiality_datafile": f"{data_path}/essential_genes.csv"
    }

    output_dir = Path(scratch)
    classifier_dir = data_path

    # Initialize KBDataLakeUtils for genome conversion and table building
    kbdl = KBDataLakeUtils(
        reference_path=sdk_config["reference_path"],
        module_path=sdk_config["module_path"],
        kb_version=kbversion,
        token=token
    )

    # Step 1: Convert user genomes to full-format TSV from workspace
    user_genome_dir = output_dir / "genome"
    user_genome_dir.mkdir(parents=True, exist_ok=True)

    for input_ref in input_refs:
        info = kbdl.get_object_info(input_ref)
        genome_name = info[1]
        user_genome_tsv = user_genome_dir / f"user_{genome_name}_kbasedump.tsv"
        print(f"Converting genome {input_ref} to TSV...")
        kbdl.run_user_genome_to_tsv(input_ref, str(user_genome_tsv))
        print(f"  Wrote: {user_genome_tsv}")

    # Step 2: Collect RAST-annotated TSVs for model reconstruction
    _skip_suffixes = ("_kbasedump.tsv", "_bakta.tsv", "_KOfamscan.tsv", "_PSORT.tsv", "_genome_data.tsv")
    all_tsvs = [f for f in sorted(user_genome_dir.glob("user_*.tsv"))
                if not f.name.endswith(_skip_suffixes)]
    print(f"Found {len(all_tsvs)} user genome TSVs in {user_genome_dir}")
    print(all_tsvs)

    # Add pangenome member TSVs
    """
    pangenome_dir = output_dir / "pangenome"
    if pangenome_dir.exists():
        for clade_dir in pangenome_dir.iterdir():
            if clade_dir.is_dir():
                genome_dir = clade_dir / "genome"
                if genome_dir.exists():
                    rast_tsvs = [f for f in sorted(genome_dir.glob("*.tsv"))
                                 if not f.name.endswith(_skip_suffixes)]
                    all_tsvs.extend(rast_tsvs)
                    print(f"Found {len(rast_tsvs)} pangenome TSVs in {genome_dir}")

    print(f"Total TSVs to process: {len(all_tsvs)}")
    print(all_tsvs)
    """
    all_tsvs = []
    pangenome_dir = output_dir / "pangenome"
    genomes_to_process = {}
    filename_prefix_rast = '_rast.tsv'  # RAST annotation prefix
    for filename_rast in user_genome_dir.glob('*' + filename_prefix_rast):
        all_tsvs.append(filename_rast)
    for filename_rast in pangenome_dir.glob('**/*' + filename_prefix_rast):
        all_tsvs.append(filename_rast)
    for filename_rast in all_tsvs:
        genome_id = filename_rast.name[:-len(filename_prefix_rast)]  # get genome_id
        print(genome_id, filename_rast)
        genomes_to_process[genome_id] = (
            read_rast_as_genome(filename_rast, genome_id),
            output_dir / 'models' / genome_id)

    # Step 3: Run model reconstruction in parallel
    results_dir = output_dir / "models"
    results_dir.mkdir(parents=True, exist_ok=True)

    work_items = []
    for tsv_path in all_tsvs:
        stem = tsv_path.stem
        output_base = str(results_dir / stem)
        work_items.append((str(tsv_path), output_base))

    print(f"\n--- Model Reconstruction ({len(work_items)} genomes, 10 workers) ---")
    futures = {}
    with ProcessPoolExecutor(max_workers=10) as executor:
        for genome_id, (genome, outp) in genomes_to_process.items():
            print(f'submit - run_model_reconstruction {genome_id} {genome} {outp} {classifier_dir} {kbversion}')
            _future = executor.submit(run_model_reconstruction2, genome_id, genome, str(outp), classifier_dir, kbversion, token)
            futures[_future] = (genome_id, outp)
        for future in as_completed(futures):
            genome_id, outp = futures[future]
            try:
                result = future.result()
                status = "OK" if result.get('success') else f"FAIL: {str(result.get('error', '?'))[:80]}"
                print(f"  {genome_id}: {status}")
            except Exception as e:
                print(f"  {genome_id}: ERROR - {e}")
                traceback.print_exc()
        """
        for inp, outp in work_items:
            print(f'submit - run_model_reconstruction {inp} {outp} {classifier_dir} {kbversion}')
            _future = executor.submit(run_model_reconstruction, inp, outp, classifier_dir, kbversion,token)
            futures[_future] = (inp, outp)

        for future in as_completed(futures):
            inp, outp = futures[future]
            try:
                result = future.result()
                status = "OK" if result.get('success') else f"FAIL: {str(result.get('error', '?'))[:80]}"
                print(f"  {Path(inp).name}: {status}")
            except Exception as e:
                print(f"  {Path(inp).name}: ERROR - {e}")
                traceback.print_exc()
        """

    cobra_files = sorted(results_dir.glob("*_cobra.json"))
    data_files = sorted(results_dir.glob("*_data.json"))
    print(f"\nModel reconstruction done. {len(cobra_files)} cobra models, {len(data_files)} data files")
    # Step 4: Run phenotype simulation in parallel
    pheno_items = []
    phenopath = output_dir / "phenotypes"
    phenopath.mkdir(parents=True, exist_ok=True)

    for cobra_path in cobra_files:
        genome_id = cobra_path.name.replace("_cobra.json", "")
        pheno_output = str(phenopath / f"{genome_id}_phenosim.json")
        pheno_items.append((str(cobra_path), pheno_output))

    print(f"\n--- Phenotype Simulation ({len(pheno_items)} models, 10 workers) ---")
    with ProcessPoolExecutor(max_workers=10) as executor:
        futures = {
            executor.submit(run_phenotype_simulation, model_file, pheno_file, data_path, max_phenotypes, kbversion,token): model_file
            for model_file, pheno_file in pheno_items
        }
        for future in as_completed(futures):
            model_path = futures[future]
            try:
                result = future.result()
                status = "OK" if result.get('success') else f"FAIL: {str(result.get('error', '?'))[:80]}"
                print(f"  {Path(model_path).name}: {status}")
            except Exception as e:
                print(f"  {Path(model_path).name}: ERROR - {e}")

    pheno_files = sorted(phenopath.glob("*_phenosim.json"))
    print(f"\nPhenotype simulation done. {len(pheno_files)} phenotype files")

    # Step 5: Build phenotype tables
    print(f"\nBuilding phenotype tables...")
    kbdl.build_phenotype_tables(
        output_dir=str(phenopath),
        phenosim_directory=str(phenopath),
        experiment_data_file=sdk_config["experimental_phenotype_datafile"],
        phenoset_file=sdk_config["phenotypeset_file"],
        fitness_mapping_dir=str(output_dir / "genome"),
        fitness_genomes_dir=sdk_config["fitness_genomes_dir"],
        model_data_dir=str(results_dir),
        reference_phenosim_dir=sdk_config["reference_phenosim_dir"],
        essential_genes_file=sdk_config["experimental_essentiality_datafile"],
    )

    # Step 6: Build model tables
    print(f"\nBuilding model tables...")
    kbdl.build_model_tables(model_path=str(results_dir))

    # Step 7: Build BERDL database - DO NOT USE THIS CODE - FILIPE IS WRITING THE CODE FOR THIS
    #db_path = output_dir / "berdl_tables.db"
    #print(f"\nBuilding BERDL database: {db_path}")
    #build_berdl_database(
    #    str(output_dir),
    #    str(db_path),
    #    reference_data_path=sdk_config["reference_path"],
    #)

    print(f"\nModel pipeline complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run model reconstruction, phenotype simulation, and table building pipeline"
    )
    parser.add_argument(
        "input_params",
        help="Path to input params JSON file"
    )
    args = parser.parse_args()

    if not os.path.exists(args.input_params):
        raise FileNotFoundError(f"Input params file not found: {args.input_params}")

    with open(args.input_params, 'r') as fh:
        _params = json.load(fh)

    main(_params)
