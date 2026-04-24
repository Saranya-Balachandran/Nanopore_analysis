// Submit example: nextflow run count_to_seurat5.nf -params-file count_to_seurat5_params.yaml
// include modules
include { basecalling } from './basecalling_modules.nf'
include { demux } from './basecalling_modules.nf'


/*--------------MAIN WORKFLOW----------------*/

workflow {

    // create a channel for parallelising the processing of all samples
    samples_channel=Channel.fromPath(params.samplesheet)
    | splitCsv( header: true )
    | map { row ->
        [row.unique_name, row.pod5_path, row.outfolder]
    }
    samples_channel.view()
    basecalling(samples_channel)
    demux(basecalling.out)

  }
