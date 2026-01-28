/*
A KBase module: KBDatalakeApps

This module provides applications for interacting with the KBase data lake,
including data retrieval, processing, and analysis utilities.

Author: chenry
*/

module KBDatalakeApps {

    /*
    Standard report output structure used by KBase apps
    */
    typedef structure {
        string report_name;
        string report_ref;
        string workspace;
    } ReportResults;

    /*
    Parameters for building genome datalake tables

    input_refs - list of workspace references to Genome or GenomeSet objects
    suffix - string suffix to append to generated table names
    save_models - boolean (0/1) indicating if generated models should be saved
    workspace_name - the name of the workspace for saving results
    */
    typedef structure {
        list<string> input_refs;
        string suffix;
        int save_models;
        string workspace_name;
    } BuildGenomeDatalakeTablesParams;

    /*
    Build genome datalake tables from Genome or GenomeSet objects

    This function takes a list of Genome or GenomeSet references and builds
    datalake tables from them. Optionally saves generated models to the workspace.
    */
    funcdef build_genome_datalake_tables(BuildGenomeDatalakeTablesParams params)
        returns (ReportResults output)
        authentication required;

};
