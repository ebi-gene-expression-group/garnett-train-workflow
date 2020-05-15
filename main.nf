#!/usr/bin/env nextflow 

// build training CDS object 
TRAINING_10X_DIR = Channel.fromPath(params.training_10x_dir)
process build_train_CDS_object{
    conda "${baseDir}/envs/monocle3-cli.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    input: 
        file(training_10x_dir) from TRAINING_10X_DIR

    output:
        file("training_cds.rds") into TRAINING_CDS
    
    """
    monocle3 create training_cds.rds\
                --expression-matrix ${training_10x_dir}/matrix.mtx\
                --cell-metadata ${training_10x_dir}/barcodes.tsv\
                --gene-annotation ${training_10x_dir}/genes.tsv\
    """
}
TRAINING_CDS = TRAINING_CDS.first()


// transform markers from SCXA format into Garnett
SCXA_MARKER_GENES = Channel.fromPath(params.marker_genes).first()
process transform_markers{
    conda "${baseDir}/envs/garnett-cli.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    input:
        file(scxa_markers) from SCXA_MARKER_GENES

    output:
        file("garnett_markers.txt") into GARNETT_MARKERS
        file("markers_list.rds") into MARKERS_LIST

    """
    transform_marker_file.R\
            --input-marker-file ${scxa_markers}\
            --pval-col ${params.pval_col}\
            --gene-names ${params.gene_names}\
            --groups-col ${params.groups_col}\
            --marker-list markers_list.rds\
            --garnett-marker-file garnett_markers.txt
    """
}

// check supplied markers 
process check_markers{ 
    publishDir "data/output_dir", mode: 'copy'
    conda "${baseDir}/envs/garnett-cli.yaml" 
    
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    input:
        file(marker_genes) from GARNETT_MARKERS
        file(training_cds) from TRAINING_CDS

    output:
        file("marker_genes_checked.txt") into CHECKED_MARKER_GENES
        

    """
    garnett_check_markers.R\
            --cds-object ${training_cds}\
            --marker-file-path ${marker_genes}\
            --marker-file-gene-id-type ${params.marker_gene_id_type}\
            -d ${params.database}\
            --cds-gene-id-type ${params.training_cds_gene_id_type}\
            --marker-output-path marker_genes_checked.txt\
    """
}

process update_markers {

    conda "${baseDir}/envs/garnett-cli.yaml" 
    
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    input:
        file(marker_list) from MARKERS_LIST
        file(marker_summary) from CHECKED_MARKER_GENES

    output: 
        file("garnett_markers_upd.txt") into GARNETT_MARKERS_UPD

    """
    update_marker_file.R\
            --marker-list-obj ${marker_list}\
            --marker-check-file ${marker_summary}\
            --updated-marker-file garnett_markers_upd.txt
    """

}

process train_classifier{
    publishDir "${params.results_dir}", mode: 'copy'
    conda "${baseDir}/envs/garnett-cli.yaml" 

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    input:
        file(training_cds) from TRAINING_CDS
        file(marker_genes) from GARNETT_MARKERS_UPD

    output:
        file("garnett_classifier.rds") into TRAINED_CLASSIFIER 

    """
    garnett_train_classifier.R\
            --cds-object ${training_cds}\
            --marker-file-path ${marker_genes}\
            -d ${params.database}\
            --train-id ${params.training_dataset_id}\
            --cds-gene-id-type ${params.training_cds_gene_id_type}\
            --marker-file-gene-id-type ${params.marker_gene_id_type}\
            --classifier-gene-id-type ${params.classifier_gene_type}\
            -n ${params.n_outgroups}\
            --output-path trained_classifier.rds 
    """
}

