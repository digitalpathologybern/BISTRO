nextflow.enable.dsl = 2

// ============================================================================
// EXISTING PROCESSES (unchanged)
// ============================================================================

// 1. Read the filtered.zarr and save a temporary .csv with the cellxgene expression
process readZarr {
    label 'python_process'

    input:
    path zarrFile
    val outputDir

    output:
    path "$outputDir/*_counts.csv"
    path "$outputDir/*_metadata.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_counts.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_metadata.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_counts.csv $outputDir/
        ln -s "\$PUB"/*_metadata.csv $outputDir/
        echo "Restored from published output"
    else
        python ${projectDir}/bin/save_expression_matrix.py $zarrFile $outputDir filtered
    fi
    """
}

// 1b. Assign FOV identifiers (centralized, runs once)
//     For CosMx: verifies native FOVs, computes FOV centers.
//     For Xenium/MERSCOPE: generates pseudo-FOVs via rasterization.
process assignFOV {
    label 'python_process'

    input:
    path metadataCSV
    val technology
    val outputDir

    output:
    path "$outputDir/*_metadata.csv"
    path "$outputDir/fov_info.json"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_metadata.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/fov_info.json 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_metadata.csv $outputDir/
        ln -s "\$PUB"/fov_info.json $outputDir/
        echo "Restored from published output"
    else
        export PYTHONPATH=${projectDir}/utils:\${PYTHONPATH:-}
        python ${projectDir}/bin/preprocessing/assign_fov.py \\
            --metadata $metadataCSV \\
            --technology $technology \\
            --output_dir $outputDir
    fi
    """
}
// 2. Apply normalization methods
//// 2.1 Apply the SCTransform to the count data
process scTransformV2 {
    label 'r_process'

    input: 
    path countMatrix
    val outputDir

    output:
    path "$outputDir/*_scTransform.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_scTransform*.csv 1>/dev/null 2>&1; then
        for f in "\$PUB"/*_scTransform*.csv; do
            ln -s "\$f" $outputDir/
        done
        echo "Restored from published output"
    else
        Rscript ${projectDir}/bin/normalization/scTransform_norm.R $countMatrix $outputDir
    fi
    """
}

//// 2.2 Apply the Area Normalization to the count data
process areaNorm {
    label 'python_process'

    input: 
    path countMatrix
    path metadataMatrix
    val outputDir

    output:
    path "$outputDir/*_areaNorm.csv"
    path "$outputDir/*_areaSF.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_areaNorm.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_areaSF.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_areaNorm.csv $outputDir/
        ln -s "\$PUB"/*_areaSF.csv $outputDir/
        echo "Restored from published output"
    else
        python ${projectDir}/bin/normalization/area_norm.py $countMatrix $metadataMatrix $outputDir
    fi
    """
}

//// 2.3 Apply cpm normalization
process cpmNorm {
    label 'r_process'

    input: 
    path countMatrix
    val outputDir

    output:
    path "$outputDir/*_cpm.csv"
    path "$outputDir/*_cpmSF.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_cpm.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_cpmSF.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_cpm.csv $outputDir/
        ln -s "\$PUB"/*_cpmSF.csv $outputDir/
        echo "Restored from published output"
    else
        Rscript ${projectDir}/bin/normalization/ls_norm.R $countMatrix $outputDir 1000000
    fi
    """
}

//// 2.4. Apply cp10k normalization
process cp10kNorm {
    label 'r_process'

    input: 
    path countMatrix
    val outputDir

    output:
    path "$outputDir/*_cp10k.csv"
    path "$outputDir/*_cp10kSF.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_cp10k.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_cp10kSF.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_cp10k.csv $outputDir/
        ln -s "\$PUB"/*_cp10kSF.csv $outputDir/
        echo "Restored from published output"
    else
        Rscript ${projectDir}/bin/normalization/ls_norm.R $countMatrix $outputDir 10000
    fi
    """
}

//// 2.5. Apply cp100 normalization
process cp100Norm {
    label 'r_process'

    input: 
    path countMatrix
    val outputDir

    output:
    path "$outputDir/*_cp100.csv"
    path "$outputDir/*_cp100SF.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_cp100.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_cp100SF.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_cp100.csv $outputDir/
        ln -s "\$PUB"/*_cp100SF.csv $outputDir/
        echo "Restored from published output"
    else
        Rscript ${projectDir}/bin/normalization/ls_norm.R $countMatrix $outputDir 100
    fi
    """
}


//// 2.6 Apply tmm normalization
process tmmNorm {
    label 'r_process'

    input: 
    path countMatrix
    val outputDir

    output:
    path "$outputDir/*_tmm.csv"
    path "$outputDir/*_tmmSF.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_tmm.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_tmmSF.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_tmm.csv $outputDir/
        ln -s "\$PUB"/*_tmmSF.csv $outputDir/
        echo "Restored from published output"
    else
        Rscript ${projectDir}/bin/normalization/tmm_norm.R $countMatrix $outputDir
    fi
    """
}

//// 2.7 Apply scran normalization
process scranNorm {
    label 'r_process'

    input: 
    path countMatrix
    val outputDir

    output:
    path "$outputDir/*_scran.csv"
    path "$outputDir/*_scranSF.csv"


    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_scran.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_scranSF.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_scran.csv $outputDir/
        ln -s "\$PUB"/*_scranSF.csv $outputDir/
        echo "Restored from published output"
    else
        Rscript ${projectDir}/bin/normalization/scran_norm.R $countMatrix $outputDir
    fi
    """
}

//// 2.8 Apply deseq2 normalization
process deseq2Norm {
    label 'r_process'

    input: 
    path countMatrix
    val outputDir

    output:
    path "$outputDir/*_deseq2.csv"
    path "$outputDir/*_deseq2SF.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_deseq2.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_deseq2SF.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_deseq2.csv $outputDir/
        ln -s "\$PUB"/*_deseq2SF.csv $outputDir/
        echo "Restored from published output"
    else
        Rscript ${projectDir}/bin/normalization/deseq2_norm.R $countMatrix $outputDir
    fi
    """

}

//// 2.9 Apply SpaNorm normalization (logpac)
process spaNormLog {
    label 'r_process'


    input:
    path countMatrix
    path metadataMatrix
    val outputDir
    val separate_fovs
    val technology

    output:
    path "$outputDir/*_spanorm-logpac.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    export SPANORM_CHECKPOINT_DIR="${params.output_folder_path}/$outputDir"
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_spanorm-logpac.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_spanorm-logpac.csv $outputDir/
        echo "Restored from published output"
    else
        Rscript ${projectDir}/bin/normalization/spanorm.R $countMatrix $metadataMatrix $outputDir $separate_fovs $technology "logpac"
    fi
    """
}

//// 2.10 Apply SpaNorm normalization (pearson)
process spaNormPearson {
    label 'r_process'

    
    input:
    path countMatrix
    path metadataMatrix
    val outputDir
    val separate_fovs
    val technology

    output:
    path "$outputDir/*_spanorm-pearson.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    export SPANORM_CHECKPOINT_DIR="${params.output_folder_path}/$outputDir"
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_spanorm-pearson.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_spanorm-pearson.csv $outputDir/
        echo "Restored from published output"
    else
        Rscript ${projectDir}/bin/normalization/spanorm.R $countMatrix $metadataMatrix $outputDir $separate_fovs $technology "pearson"
    fi
    """
}

//// 2.11 Apply no normalization
process noneNorm {
    label 'r_process'

    input:
    path countMatrix
    val outputDir

    output:
    path "$outputDir/*_none.csv"
    path "$outputDir/*_noneSF.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_none.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_noneSF.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*_none.csv $outputDir/
        ln -s "\$PUB"/*_noneSF.csv $outputDir/
        echo "Restored from published output"
    else
        Rscript ${projectDir}/bin/normalization/none_norm.R $countMatrix $outputDir
    fi
    """
}

// 4. Select the HVG for each normalization method
process hvg_selection {
    label 'python_process'

    input:
    path norm_csv
    val outputHVG

    output:
    path "$outputHVG/*.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputHVG
    python ${projectDir}/bin/find_hvg.py --norm_counts_csv $norm_csv --output_csv $outputHVG
    """
}

// 5. Load all HVG CSV files and generate upset plots
process hvg_upset_plot {
    label 'python_process'

    input:
    path hvgCSVList
    val outputPlot
    val outputCSV

    output:
    path "$outputPlot/*.png"
    path "$outputCSV/top_*.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputPlot
    mkdir -p $outputCSV
    PUB_PLOT="${params.output_folder_path}/$outputPlot"
    PUB_CSV="${params.output_folder_path}/$outputCSV"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB_PLOT"/*.png 1>/dev/null 2>&1 && \
       ls "\$PUB_CSV"/*.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB_PLOT"/*.png $outputPlot/
        ln -s "\$PUB_CSV"/top_*.csv $outputCSV/
        echo "Restored from published output"
    else
        python ${projectDir}/bin/plotters/hvg_upset_plot.py --hvg_csv $hvgCSVList --top_fraction ${params.hvg_fractions} --output_plot $outputPlot --output_csv $outputCSV
    fi
    """
}

// 6. Create a pathway analysis based on the HVG
process pathway_analysis_hvg {
    label 'r_process'

    input:
    path topHVGCSV
    val outputPlot
    val outputPathway

    output:
    path "$outputPlot/pathway/*.png"
    path "$outputPathway/*.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    module load CMake
    mkdir -p $outputPlot
    mkdir -p $outputPlot/pathway
    mkdir -p $outputPathway
    Rscript ${projectDir}/bin/plotters/pathway_analysis_HVG.R $topHVGCSV $outputPlot/pathway $outputPathway
    """
}

// 7. Generate a GSVA analysis based on the HVG for each norm method
process gsva_analysis {
    label 'r_process'

    input:
    path topHVGCSV
    path expressionMatricesCSV
    val outputPlot

    output:
    path "$outputPlot/gsva/*.png"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputPlot
    mkdir -p $outputPlot/gsva
    Rscript ${projectDir}/bin/plotters/gsva_plot.R $topHVGCSV $outputPlot/gsva $expressionMatricesCSV
    """
}

// 8. Run single cell annotation using InSituType
process insitutype_annotation {
    label 'r_process'

    input:
    path referenceAtlas
    path countMatrix
    path normMatrices
    path metadataMatrix
    path hvgCSV
    val outputDir

    output:
    path "$outputDir/*.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    FRAC=\$(echo $hvgCSV | grep -oP 'top_\\K[0-9.]+')
    if [ "${params.restore_published}" = "true" ] && \
       [ -n "\$FRAC" ] && \
       ls "\$PUB"/\${FRAC}_HVG_*_annotation.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/\${FRAC}_HVG_*_annotation.csv $outputDir/
        echo "Restored fraction \$FRAC from published output"
    else
        Rscript ${projectDir}/bin/annotation/run_insitutype.R $referenceAtlas $countMatrix $metadataMatrix $hvgCSV $outputDir $normMatrices
    fi
    """
}

process cell_type_to_colormap {
    label 'r_process'

    input:
    path cellTypeCSV
    val color_list

    output:
    path "color_palette.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    Rscript ${projectDir}/bin/plotters/create_colormap.R $color_list $cellTypeCSV
    """
}


// ============================================================================
// NEW EVALUATION PROCESSES
// ============================================================================

// 10. FOV Batch Effect Evaluation (Section 2.2)
//     Runs AFTER all normalizations are complete.
//     Fits OLS and MELM per normalization, outputs summary CSV,
//     per-FOV random intercepts, and serialized model objects.
process batch_effect_evaluation {
    label 'python_process'

    input:
    path zarrFile
    val nextflowOutputPath
    val technology
    val tissueAnnotation
    val outputDir
    val datasetName
    path enrichedMetadata
    val heAlignmentPath
    val pixelSize
    path _normDone

    output:
    path "$outputDir/*_batch_effect_summary.csv"
    path "$outputDir/*_random_intercepts.csv"
    path "$outputDir/*_batch_effect_models.pkl"
    path "$outputDir/*_fov_summary.csv"
    path "$outputDir/*_tissue_ls.csv", optional: true
    path "$outputDir/*_drift_comparison.csv", optional: true
    path "$outputDir/*_drift_intercepts.csv", optional: true
    path "$outputDir/*_fov_centers_*.csv", optional: true
    path "$outputDir/*_cell_overview.csv", optional: true

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    def he_arg = heAlignmentPath ? "--he_alignment \"${heAlignmentPath}\"" : ""
    def vc_arg = params.vc_column ? "--vc_column ${params.vc_column}" : ""
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_batch_effect_summary.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_random_intercepts.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_batch_effect_models.pkl 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_fov_summary.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*.csv $outputDir/ 2>/dev/null
        ln -s "\$PUB"/*.pkl $outputDir/ 2>/dev/null
        echo "Restored from published output"
    else
        export PYTHONPATH=${projectDir}/utils:\${PYTHONPATH:-}
        python ${projectDir}/bin/evaluation/batch_effect.py \\
            --zarr $zarrFile \\
            --nextflow_output "$nextflowOutputPath" \\
            --technology $technology \\
            --tissue_annotation "$tissueAnnotation" \\
            --output_dir $outputDir \\
            --dataset_name "$datasetName" \\
            --metadata $enrichedMetadata \\
            --pixel_size $pixelSize \\
            --n_boot_ci ${params.n_boot_ci} \\
            $he_arg \\
            $vc_arg \\
            --use_log
    fi
    """
}

// 11. Transformation / Variance Stabilization Analysis (Section 2.3)
//     Runs AFTER all normalizations are complete.
//     Applies transformations, computes mean-variance, PCA, binned CV,
//     PC1 correlation for each normalization-transformation combination.
process transformation_analysis {
    label 'python_process'

    input:
    path zarrFile
    val nextflowOutputPath
    val outputDir
    val datasetName
    val technology
    path _normDone       // dummy input to enforce ordering after normalization

    output:
    path "$outputDir/*_mean_variance_stats.csv"
    path "$outputDir/*_transformation_summary.csv"
    path "$outputDir/*_pc1_correlation.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_mean_variance_stats.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_transformation_summary.csv 1>/dev/null 2>&1 && \
       ls "\$PUB"/*_pc1_correlation.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*.csv $outputDir/
        echo "Restored from published output"
    else
        export PYTHONPATH=${projectDir}/utils:\${PYTHONPATH:-}
        python ${projectDir}/bin/evaluation/transformation_analysis.py \\
            --zarr "$zarrFile" \\
            --nextflow_output "$nextflowOutputPath" \\
            --output_dir "$outputDir" \\
            --dataset_name "$datasetName" \\
            --technology "$technology" \\
            --pseudocounts 0.01 0.1 0.5 1 10 \\
            --alpha 0.05 \\
            --theta 100
    fi
    """
}

// 12. HVG Selection Benchmark (Section 2.4)
//     Runs PER LAYER after HVG selection is complete.
//     Evaluates clustering stability, HVG vs random, spatial coherence, AUROC.
//     Parallelized across layers by Nextflow.
process hvg_benchmark {
    label 'python_process'

    input:
    path zarrFile
    val nextflowOutputPath
    val technology
    val datasetName
    val layerName
    val hvgFraction
    val outputDir
    path referenceAnnotation
    path _hvgDone

    output:
    path "$outputDir/checkpoints/*.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir/checkpoints
    PUB="${params.output_folder_path}/$outputDir"
    FRAC=\$(printf "%.2f" ${hvgFraction})
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/checkpoints/${layerName}_\${FRAC}.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/checkpoints/${layerName}_\${FRAC}.csv $outputDir/checkpoints/
        echo "Restored ${layerName}_\${FRAC} from published output"
    else
        export PYTHONPATH=${projectDir}/utils:\${PYTHONPATH:-}
        python ${projectDir}/bin/evaluation/hvg_benchmark.py \\
            --zarr $zarrFile \\
            --nextflow_output "$nextflowOutputPath" \\
            --technology $technology \\
            --dataset_name "$datasetName" \\
            --layer $layerName \\
            --output_dir $outputDir \\
            --checkpoint_dir "${params.output_folder_path}/${params.outputEval}/hvg_benchmark" \\
            --n_bootstrap ${params.hvg_n_bootstrap} \\
            --subsample_frac ${params.hvg_subsample_frac} \\
            --n_random_repeats ${params.hvg_n_random_repeats} \\
            --k_neighbors ${params.hvg_k_neighbors} \\
            --resolution ${params.hvg_resolution} \\
            --seed ${params.hvg_seed} \\
            --hvg_fractions $hvgFraction \\
            --reference_annotation $referenceAnnotation \\
            --n_bootstrap_phase3 ${params.hvg_n_bootstrap}
    fi
    """
}

// 13. Cell Type Annotation Agreement (Section 2.4, Figure 4D)
//     Runs AFTER InSituType annotation is complete.
//     Computes pairwise ARI between normalization methods.
process annotation_agreement {
    label 'python_process'

    input:
    val annotationDir
    val hvgThreshold
    val outputDir
    val datasetName
    path _annoDone       // dummy input to enforce ordering after annotation

    output:
    path "$outputDir/*_pairwise_ari_*.csv"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    mkdir -p $outputDir
    PUB="${params.output_folder_path}/$outputDir"
    if [ "${params.restore_published}" = "true" ] && \
       ls "\$PUB"/*_pairwise_ari_*.csv 1>/dev/null 2>&1; then
        ln -s "\$PUB"/*.csv $outputDir/
        echo "Restored from published output"
    else
        export PYTHONPATH=${projectDir}/utils:\${PYTHONPATH:-}
        python ${projectDir}/bin/evaluation/annotation_agreement.py \\
            --annotation_dir "$annotationDir" \\
            --hvg_threshold $hvgThreshold \\
            --output_dir "$outputDir" \\
            --dataset_name $datasetName
    fi
    """
}

// 14. Generate HTML Report
//     Runs LAST, after all evaluation steps are complete.
//     Collects all CSVs and produces a self-contained HTML report.
process generate_report {
    label 'python_process'

    input:
    val resultsDir
    val datasetName
    val outputHtml
    path _batchDone
    path _transDone
    path _hvgBenchDone
    path _ariDone

    output:
    path "$outputHtml"

    publishDir "${params.output_folder_path}", mode: 'copy'

    script:
    """
    export PYTHONPATH=${projectDir}/utils:\${PYTHONPATH:-}
    python ${projectDir}/bin/report/generate_report.py \\
        --results_dir "$resultsDir" \\
        --dataset_name "$datasetName" \\
        --output_html $outputHtml
    """
}


// ============================================================================
// WORKFLOW
// ============================================================================
params.restore_published = params.restore_published ?: true

workflow {

    // ---- Parameters (existing) ----
    outputNorm      = params.outputNorm
    outputPlots     = params.outputPlots
    outputTrans     = params.outputTrans
    outputHVG       = params.outputHVG
    outputDims      = params.outputDims
    outputAnno      = params.outputAnno
    outputPathway   = params.outputPathway
    temporary       = params.temporary
    color_list      = params.color_list
    zarrFile        = file(params.zarrFile)
    separate_fovs   = params.separate_fovs
    add_global      = params.add_global
    technology      = params.technology
    heAlignmentPath = params.heAlignmentPath
    pixelSize       = params.pixelSize
    // rasterize_image is removed (handled by assignFOV)

    // ---- Parameters (evaluation) ----
    outputEval          = params.outputEval
    tissueAnnotation    = params.tissueAnnotation
    datasetName         = params.datasetName
    nextflowOutputPath  = params.output_folder_path
    hvgThreshold        = params.hvgThreshold ?: '0.75'

    // ---- Toggle flags ----
    skip_annotation = params.skip_annotation ?: false
    skip_pathway    = params.skip_pathway ?: false
    skip_hvg_bench    = params.skip_hvg_bench ?: false

    // ========================================================================
    // STEP 1: Read zarr
    // ========================================================================
    (countMatrix, rawMetadata) = readZarr(zarrFile, outputNorm)

    // ========================================================================
    // STEP 1b: Assign FOV identifiers (centralized)
    // ========================================================================
    (metadataMatrix, fovInfo) = assignFOV(rawMetadata, technology, outputNorm)

    // ========================================================================
    // STEP 2: Normalization (all methods in parallel)
    // ========================================================================
    sctransformFile = scTransformV2(countMatrix, outputNorm)
    (areaNormFile, areaSFFile) = areaNorm(countMatrix, metadataMatrix, outputNorm)
    (cpmNormFile, cpmSFFile) = cpmNorm(countMatrix, outputNorm)
    (cp10kNormFile, cp10kSFFile) = cp10kNorm(countMatrix, outputNorm)
    (cp100NormFile, cp100SFFile) = cp100Norm(countMatrix, outputNorm)
    (tmmNormFile, tmmSFFile) = tmmNorm(countMatrix, outputNorm)
    (scranNormFile, scranSFFile) = scranNorm(countMatrix, outputNorm)
    (deseq2NormFile, deseq2SFFile) = deseq2Norm(countMatrix, outputNorm)
    spaNormFileLog = spaNormLog(countMatrix, metadataMatrix, outputNorm, separate_fovs, technology)
    spaNormFilePears = spaNormPearson(countMatrix, metadataMatrix, outputNorm, separate_fovs, technology)
    (noneFile, noneSFFile) = noneNorm(countMatrix, outputNorm)

    // ========================================================================
    // STEP 3: Collect normalization outputs
    // ========================================================================
    sfnormalizedCSVChannel = areaNormFile.concat(cpmNormFile, cp10kNormFile, cp100NormFile, tmmNormFile, scranNormFile, deseq2NormFile, noneFile)
    allNormDataCSVChannel = sfnormalizedCSVChannel.concat(sctransformFile, spaNormFileLog, spaNormFilePears)
    allNormCollected = allNormDataCSVChannel.collect()

    // ========================================================================
    // STEP 4: HVG selection
    // ========================================================================
    hvgCSVChannel = hvg_selection(allNormDataCSVChannel, outputHVG)

    // ========================================================================
    // STEP 5: HVG upset plots
    // ========================================================================
    (_, selectedHVGCSV) = hvg_upset_plot(hvgCSVChannel.collect(), outputPlots, outputHVG)
    hvgCollected = selectedHVGCSV.collect()

    // ========================================================================
    // STEP 6: Pathway analysis
    // ========================================================================
    if (!skip_pathway) {
        (_, _) = pathway_analysis_hvg(selectedHVGCSV.flatMap(), outputPlots, outputPathway)
    }

    // ========================================================================
    // STEP 7: GSVA (uncomment to enable)
    // ========================================================================
    //_ = gsva_analysis(selectedHVGCSV.flatMap(), allNormDataCSVChannel.collect(), outputPlots)

    // ========================================================================
    // STEP 8: Cell type annotation (InSituType)
    // ========================================================================
    if (!skip_annotation) {
        scReferenceFile = file(params.scReferenceFile)
        annotations = insitutype_annotation(scReferenceFile, countMatrix, allNormCollected, metadataMatrix, selectedHVGCSV.flatMap(), outputAnno)
        palette_csv = cell_type_to_colormap(annotations.collect(), color_list)
        annoCollected = annotations.collect()

        // Reference annotation for HVG benchmark (from InSituType "none" at 100% HVG)
        referenceAnnotationPath = annotations.flatten().filter { it.name ==~ /.*1.0_HVG_none_annotation\.csv/ }.first()
    }


    // ========================================================================
    // STEP 10: FOV Batch Effect Evaluation
    // ========================================================================
    (batchSummary, batchIntercepts, batchModels, batchFovSummary, batchTissueLs, batchDriftComp, batchDriftInt, batchDriftCenters, batchCellOverview) = batch_effect_evaluation(
        zarrFile,
        nextflowOutputPath,
        technology,
        tissueAnnotation,
        outputEval + '/batch_effect',
        datasetName,
        metadataMatrix,
        heAlignmentPath,
        pixelSize,
        allNormCollected
    )

    // ========================================================================
    // STEP 11: Transformation / Variance Stabilization Analysis
    // ========================================================================
    (transMeanVar, transSummary, transPC1) = transformation_analysis(
        zarrFile,
        nextflowOutputPath,
        outputEval + '/transformation',
        datasetName,
        technology,
        allNormCollected
    )

    // ========================================================================
    // STEP 12: HVG Selection Benchmark (per layer x fraction, parallelized)
    // ========================================================================
    layerNamesChannel = allNormDataCSVChannel.map { file ->
        def basename = file.getName().replaceAll('\\.csv$', '')
        def parts = basename.split('_')
        return parts[-1]
    }

    // Create all (layer, fraction) combinations for maximum parallelism
    hvgFractionsChannel = Channel.of(params.hvg_fractions.trim().split(/\s+/)).flatten()
    layerFractionPairs = layerNamesChannel.combine(hvgFractionsChannel)

    // Reference annotation: InSituType from "none" normalization at 100% HVG
    // Extract the reference annotation file from InSituType outputs
    // (1.0_HVG_none_annotation.csv = "none" normalization at 100% HVG)
    // referenceAnnotationPath = annotations.flatten().filter { it.name ==~ /.*1.0_HVG_none_annotation\.csv/ }.first()

    // Default: empty channel when HVG benchmark is skipped, so generate_report
    // still has something to consume for its `_hvgBenchDone` input.
    hvgBenchCollected = Channel.empty().collect()

    if (!skip_hvg_bench) {
        hvgBenchResults = hvg_benchmark(
            zarrFile,
            nextflowOutputPath,
            technology,
            datasetName,
            layerFractionPairs.map { it[0] },
            layerFractionPairs.map { it[1] },
            outputEval + '/hvg_benchmark',
            referenceAnnotationPath,
            hvgCollected
        )
        hvgBenchCollected = hvgBenchResults.collect()
    }
    
        
    if (!skip_annotation) {
        // ========================================================================
        // STEP 13: Annotation Agreement
        // ========================================================================
        ariResults = annotation_agreement(
            nextflowOutputPath + '/' + outputAnno,
            hvgThreshold,
            outputEval + '/annotation_agreement',
            datasetName,
            annoCollected
        )
    

        // ========================================================================
        // STEP 14: Generate HTML Report
        // ========================================================================
        generate_report(
            nextflowOutputPath + '/' + outputEval,
            datasetName,
            'BISTRO_report.html',
            batchSummary.mix(batchIntercepts).mix(batchFovSummary).mix(batchTissueLs).mix(batchDriftComp).mix(batchDriftInt).mix(batchDriftCenters).collect(),
            transSummary.mix(transPC1).collect(),
            hvgBenchCollected,
            ariResults
        )
    }
}