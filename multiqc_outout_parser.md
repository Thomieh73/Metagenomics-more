[Multiqc](https://seqera.io/multiqc/) is a very powerful tool for checking the sequencing quality of your metagenomes or the output of many others tools that analyzed your sequencing data.

However, sometimes you want to show the results in from the fastqc process in a different way. For instance, you would like to vizualize the sequencing quality of the samples from different treatments.
That is why I created a R markdown script to process the json data file coming out of multiqc. Here is the file: [Multiqc_json_parser.Rmd](https://github.com/Thomieh73/Metagenomics-more/blob/master/Multiqc_json_parser.Rmd)

The multiqc json file was produced with the [nf-core MAG pipeline version 5.0.0](https://nf-co.re/mag/5.0.0/).
The code is a unpolished and a bit crude, but feel free to download, and modify it to your needs. And having a metadata file in tsv format with your samples helps to get your results in the way you like it. 
