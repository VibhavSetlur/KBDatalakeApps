## Output directory structure

```text
/
├── input_params.json
├── output_kbase_object.json
│
├── genome/
│   ├── user_<genome1>.faa
│   ├── user_<genome1>_rast.tsv
│   ├── user_<genome1>_kofam.tsv
│   ├── user_<genome1>_bakta.tsv
│   ├── user_<genome1>_psort.tsv
│   └── user_<genome2>.faa
│
├── assembly/
│   ├── user_<genome1>.fna
│   └── user_<genome2>.fna
│
├── library/
│   └── user_genome.txt
│
├── ani/
│   ├── kepangenomes_fast.out
│   ├── fitness_fast.out
│   └── phenotypes_fast.out
│
├── model/
│   ├── user_<genome1>/
│   └── user_<genome2>/
│
├── phenotype/
│   └── ???
│
└── pangenome/
    └── <clade_id>/
        ├── genome/
        │   ├── <member1>.faa
        │   └── <member2>.faa
        │
        ├── assembly/
        │   ├── <member1>.fna
        │   └── <member2>.fna
        │
        ├── model/
        │   ├── <member1>/
        │   └── <member2>/
        │
        ├── phenotype/
        │   └── ???
        │
        ├── mmseqs2/
        │
        ├── library_members.txt
        ├── master_faa_pangenome_members.faa
        ├── master_faa.faa
        ├── 
        └── db.sqlite
