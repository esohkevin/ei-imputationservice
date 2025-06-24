#!/usr/bin/env nextflow

def jobCompletionMessage() {
   return msg = """
          Pipeline execution summary
          ---------------------------
          Completed at: ${workflow.complete}
          Duration    : ${workflow.duration}
          Success     : ${workflow.success}
          workDir     : ${workflow.workDir}
          exit status : ${workflow.exitStatus}
          """
          .stripIndent()
}
