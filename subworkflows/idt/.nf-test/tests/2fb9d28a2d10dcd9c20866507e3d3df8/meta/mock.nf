import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test workflow
include { PROCESSFASTQ } from '/mnt/archive/work/slowney/IDT-NextflowComponents-PoC/subworkflows/idt/processfastq/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()

workflow {

    // run dependencies
    

    // workflow mapping
    def input = []
    
                skip_fastqc        = false 
                skip_downsample    = false
                umi_flag           = false
                skip_extractumi    = true
                read_structure     = '8M143T 8M143T'
                skip_trimming      = false

                // TODO nf-core: define inputs of the workflow here. Example:
                input[0] = Channel.of([
                    [ id:'test', single_end:false, sample_size:10000 ], // meta map
                    [
                        file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                        file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                    ]
                ])
                input[1] = skip_fastqc
                input[2] = skip_downsample
                input[3] = umi_flag
                input[4] = skip_extractumi
                input[5] = read_structure
                input[6] = skip_trimming
                
    //----

    //run workflow
    PROCESSFASTQ(*input)
    
    if (PROCESSFASTQ.output){

        // consumes all named output channels and stores items in a json file
        for (def name in PROCESSFASTQ.out.getNames()) {
            serializeChannel(name, PROCESSFASTQ.out.getProperty(name), jsonOutput)
        }	  
    
        // consumes all unnamed output channels and stores items in a json file
        def array = PROCESSFASTQ.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
}


def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
