"""
Ontology Term Enrichment from Local Parquet Files.

This module provides functions to enrich ontology terms (GO, EC, KEGG, COG, PFAM, SO)
with human-readable labels and definitions from local parquet files.

This is the LOCAL version that reads from files instead of calling APIs.
Use OntologyEnrichment for the API-based version.

Author: Jose P. Faria (jplfaria@gmail.com)
Date: February 2026

=== SUPPORTED ONTOLOGIES ===

- GO (Gene Ontology): GO:XXXXXXX
- EC (Enzyme Commission): EC:X.X.X.X
- KEGG (KO orthologs): KEGG:KXXXXX
- COG (Clusters of Orthologous Groups): COG:COGXXXX, COG:X (categories)
- PFAM (Protein Families): PFAM:PFXXXXX.XX
- SO (Sequence Ontology): SO:XXXXXXX

=== DATA FILES ===

Expected files at base_path (default: /data/reference_data/):
- statements.parquet: BERDL ontology statements (GO, EC, SO, PFAM)
- kegg_ko_definitions.parquet: KEGG KO definitions (ko_id, definition)
- cog_definitions.parquet: COG definitions (cog_id, category, name, gene, pathway)

=== USAGE ===

    from berdl.ontology_enrichment_local import OntologyEnrichmentLocal

    enricher = OntologyEnrichmentLocal(base_path="/data/reference_data")

    # Enrich a list of term IDs
    terms = ["GO:0008150", "EC:1.1.1.1", "KEGG:K00001"]
    enriched = enricher.enrich_terms(terms)
    # Returns: DataFrame with columns [identifier, label, definition]
"""

import pandas as pd
import pyarrow.parquet as pq
import re
from pathlib import Path
from typing import List, Dict, Optional, Set


class OntologyEnrichmentLocal:
    """
    Enrich ontology terms with labels and definitions from local parquet files.

    Data sources (local files):
    - statements.parquet: GO, EC, SO, PFAM (from BERDL ontology source)
    - kegg_ko_definitions.parquet: KEGG KO definitions
    - cog_definitions.parquet: COG definitions
    """

    # Default paths
    DEFAULT_BASE_PATH = "/data/reference_data"

    # Cache for loaded data
    _statements_df: pd.DataFrame = None
    _kegg_df: pd.DataFrame = None
    _cog_df: pd.DataFrame = None

    def __init__(self, base_path: str = None):
        """
        Initialize the local ontology enrichment client.

        Args:
            base_path: Path to directory containing parquet files.
                       Default: /data/reference_data
        """
        self.base_path = Path(base_path or self.DEFAULT_BASE_PATH)

        # File paths
        self.statements_path = self.base_path / "statements.parquet"
        self.kegg_path = self.base_path / "kegg_ko_definitions.parquet"
        self.cog_path = self.base_path / "cog_definitions.parquet"

    def _load_statements(self) -> pd.DataFrame:
        """Load statements parquet file (cached)."""
        if OntologyEnrichmentLocal._statements_df is None:
            if not self.statements_path.exists():
                raise FileNotFoundError(f"Statements file not found: {self.statements_path}")

            print(f"    Loading statements from {self.statements_path}...")
            OntologyEnrichmentLocal._statements_df = pq.read_table(self.statements_path).to_pandas()
            print(f"    Loaded {len(OntologyEnrichmentLocal._statements_df)} statements")

        return OntologyEnrichmentLocal._statements_df

    def _load_kegg(self) -> pd.DataFrame:
        """Load KEGG KO definitions parquet file (cached)."""
        if OntologyEnrichmentLocal._kegg_df is None:
            if not self.kegg_path.exists():
                raise FileNotFoundError(f"KEGG file not found: {self.kegg_path}")

            print(f"    Loading KEGG KO definitions from {self.kegg_path}...")
            OntologyEnrichmentLocal._kegg_df = pq.read_table(self.kegg_path).to_pandas()
            print(f"    Loaded {len(OntologyEnrichmentLocal._kegg_df)} KEGG KO definitions")

        return OntologyEnrichmentLocal._kegg_df

    def _load_cog(self) -> pd.DataFrame:
        """Load COG definitions parquet file (cached)."""
        if OntologyEnrichmentLocal._cog_df is None:
            if not self.cog_path.exists():
                raise FileNotFoundError(f"COG file not found: {self.cog_path}")

            print(f"    Loading COG definitions from {self.cog_path}...")
            OntologyEnrichmentLocal._cog_df = pq.read_table(self.cog_path).to_pandas()
            print(f"    Loaded {len(OntologyEnrichmentLocal._cog_df)} COG definitions")

        return OntologyEnrichmentLocal._cog_df

    def _enrich_from_statements(self, identifiers: List[str]) -> Dict[str, Dict]:
        """
        Enrich ontology terms from local statements parquet file.
        Works for GO, EC, SO, PFAM.
        """
        if not identifiers:
            return {}

        statements_df = self._load_statements()

        # Filter to relevant subjects and predicates
        mask = (
            statements_df['subject'].isin(identifiers) &
            statements_df['predicate'].isin(['rdfs:label', 'IAO:0000115'])
        )
        filtered = statements_df[mask]

        results = {}
        for _, row in filtered.iterrows():
            subj = row['subject']
            pred = row['predicate']
            val = row['value']

            if subj not in results:
                results[subj] = {'label': '', 'definition': ''}

            if pred == 'rdfs:label':
                results[subj]['label'] = val
            elif pred == 'IAO:0000115':
                results[subj]['definition'] = val

        return results

    def _enrich_kegg(self, ko_ids: List[str]) -> Dict[str, Dict]:
        """
        Enrich KEGG KO terms from local parquet file.
        """
        if not ko_ids:
            return {}

        kegg_df = self._load_kegg()

        # Create lookup dict
        ko_lookup = dict(zip(kegg_df['ko_id'], kegg_df['definition']))

        results = {}
        for ko_id in ko_ids:
            # Extract K number (e.g., K00001 from KEGG:K00001)
            k_num = ko_id.replace('KEGG:', '')

            definition = ko_lookup.get(k_num, '')
            if definition:
                # Parse label from definition
                # Format: "E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]"
                label = re.sub(r'\s*\[EC:[^\]]+\]', '', definition).strip()
                results[ko_id] = {'label': label, 'definition': definition}
            else:
                results[ko_id] = {'label': '', 'definition': ''}

        return results

    def _enrich_cog(self, cog_ids: List[str]) -> Dict[str, Dict]:
        """
        Enrich COG terms from local parquet file.
        """
        if not cog_ids:
            return {}

        cog_df = self._load_cog()

        # Create lookup dict
        cog_lookup = {}
        for _, row in cog_df.iterrows():
            cog_lookup[row['cog_id']] = {
                'name': row.get('name', ''),
                'category': row.get('category', ''),
                'gene': row.get('gene', ''),
                'pathway': row.get('pathway', '')
            }

        # COG category descriptions
        cog_categories = {
            'J': 'Translation, ribosomal structure and biogenesis',
            'A': 'RNA processing and modification',
            'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'B': 'Chromatin structure and dynamics',
            'D': 'Cell cycle control, cell division, chromosome partitioning',
            'Y': 'Nuclear structure',
            'V': 'Defense mechanisms',
            'T': 'Signal transduction mechanisms',
            'M': 'Cell wall/membrane/envelope biogenesis',
            'N': 'Cell motility',
            'Z': 'Cytoskeleton',
            'W': 'Extracellular structures',
            'U': 'Intracellular trafficking, secretion, and vesicular transport',
            'O': 'Posttranslational modification, protein turnover, chaperones',
            'X': 'Mobilome: prophages, transposons',
            'C': 'Energy production and conversion',
            'G': 'Carbohydrate transport and metabolism',
            'E': 'Amino acid transport and metabolism',
            'F': 'Nucleotide transport and metabolism',
            'H': 'Coenzyme transport and metabolism',
            'I': 'Lipid transport and metabolism',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
            'R': 'General function prediction only',
            'S': 'Function unknown',
        }

        results = {}
        for cog_id in cog_ids:
            # Handle COG categories (single letter, e.g., COG:J)
            if re.match(r'COG:[A-Z]$', cog_id):
                cat = cog_id.replace('COG:', '')
                results[cog_id] = {
                    'label': cog_categories.get(cat, ''),
                    'definition': f"COG functional category {cat}"
                }
            # Handle COG IDs (e.g., COG:COG0001)
            else:
                raw_id = cog_id.replace('COG:', '')
                info = cog_lookup.get(raw_id, {})

                if info:
                    def_parts = []
                    if info.get('category'):
                        def_parts.append(f"Category: {info['category']}")
                    if info.get('gene'):
                        def_parts.append(f"Gene: {info['gene']}")
                    if info.get('pathway'):
                        def_parts.append(f"Pathway: {info['pathway']}")

                    results[cog_id] = {
                        'label': info.get('name', ''),
                        'definition': '. '.join(def_parts) if def_parts else ''
                    }
                else:
                    results[cog_id] = {'label': '', 'definition': ''}

        return results

    def extract_ontology_terms(self, genome_df: pd.DataFrame) -> Dict[str, Set[str]]:
        """
        Extract ontology term IDs from a genome features DataFrame.

        Looks for columns like 'Annotation:GO', 'Annotation:KO', 'Annotation:EC', etc.

        Returns:
            Dict mapping ontology type to set of term IDs
        """
        terms_by_type = {
            'GO': set(),
            'EC': set(),
            'KEGG': set(),
            'COG': set(),
            'PFAM': set(),
            'SO': set(),
        }

        # Patterns for extracting term IDs
        patterns = {
            'GO': re.compile(r'GO:\d+'),
            'EC': re.compile(r'EC:[\d\.-]+'),
            'KEGG': re.compile(r'(?:KEGG:)?K\d{5}'),
            'COG': re.compile(r'COG:(?:COG\d+|[A-Z])'),
            'PFAM': re.compile(r'(?:PFAM:)?PF\d+(?:\.\d+)?'),
            'SO': re.compile(r'SO:\d+'),
        }

        # Look for annotation columns
        for col in genome_df.columns:
            if not col.startswith('Annotation:'):
                continue

            for _, row in genome_df.iterrows():
                value = str(row.get(col, ''))
                if not value or value == 'nan':
                    continue

                for ont_type, pattern in patterns.items():
                    matches = pattern.findall(value)
                    for match in matches:
                        if ont_type == 'KEGG' and not match.startswith('KEGG:'):
                            match = f'KEGG:{match}'
                        elif ont_type == 'PFAM' and not match.startswith('PFAM:'):
                            match = f'PFAM:{match}'
                        terms_by_type[ont_type].add(match)

        return terms_by_type

    def enrich_terms(self, term_ids: List[str]) -> pd.DataFrame:
        """
        Enrich a list of ontology term IDs with labels and definitions.

        Args:
            term_ids: List of ontology term IDs (e.g., ['GO:0008150', 'EC:1.1.1.1'])

        Returns:
            DataFrame with columns [identifier, label, definition]
        """
        if not term_ids:
            return pd.DataFrame(columns=['identifier', 'label', 'definition'])

        # Group by ontology type
        go_terms = [t for t in term_ids if t.startswith('GO:')]
        ec_terms = [t for t in term_ids if t.startswith('EC:')]
        kegg_terms = [t for t in term_ids if t.startswith('KEGG:')]
        cog_terms = [t for t in term_ids if t.startswith('COG:')]
        pfam_terms = [t for t in term_ids if t.startswith('PFAM:')]
        so_terms = [t for t in term_ids if t.startswith('SO:')]

        all_results = {}

        # Enrich from local statements (GO, EC, SO, PFAM)
        berdl_terms = go_terms + ec_terms + so_terms + pfam_terms
        if berdl_terms:
            print(f"  Enriching {len(berdl_terms)} terms from local statements...")
            berdl_results = self._enrich_from_statements(berdl_terms)
            all_results.update(berdl_results)

        # Enrich KEGG from local file
        if kegg_terms:
            print(f"  Enriching {len(kegg_terms)} KEGG terms from local file...")
            kegg_results = self._enrich_kegg(kegg_terms)
            all_results.update(kegg_results)

        # Enrich COG from local file
        if cog_terms:
            print(f"  Enriching {len(cog_terms)} COG terms from local file...")
            cog_results = self._enrich_cog(cog_terms)
            all_results.update(cog_results)

        # Build result DataFrame
        rows = []
        for term_id in term_ids:
            info = all_results.get(term_id, {'label': '', 'definition': ''})
            rows.append({
                'identifier': term_id,
                'label': info.get('label', ''),
                'definition': info.get('definition', '')
            })

        return pd.DataFrame(rows)

    def enrich_genome_terms(self, genome_df: pd.DataFrame) -> pd.DataFrame:
        """
        Extract and enrich all ontology terms from a genome features DataFrame.

        Args:
            genome_df: DataFrame with annotation columns (e.g., 'Annotation:GO')

        Returns:
            DataFrame with enriched ontology terms
        """
        terms_by_type = self.extract_ontology_terms(genome_df)

        all_terms = []
        for terms in terms_by_type.values():
            all_terms.extend(terms)

        if not all_terms:
            print("  No ontology terms found in genome")
            return pd.DataFrame(columns=['identifier', 'label', 'definition'])

        print(f"  Found {len(all_terms)} unique ontology terms")
        return self.enrich_terms(all_terms)

    @classmethod
    def clear_cache(cls):
        """Clear all cached data to free memory."""
        cls._statements_df = None
        cls._kegg_df = None
        cls._cog_df = None
