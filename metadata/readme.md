# readme for metadata folder

Files in this folder:

`permafrost_thaw_metadata.csv` - Not sure where this originally comes from. My notes say it was downloaded from NCBI, but I cannot find the download link again if that is where it came from. See additional data download information in folder `3 Alaska permafrost thaw shotgun sequence dataset`

`readme.md` - is this file

`samples.tsv` was created manually from `permafrost_thaw_metadata.csv`
 
 - co_assembly column allows samples to be grouped into co-assemblies for megahit
 - samples with "." should be ignored either because of sample quality or because they are control samples etc

 - note read1 and read2 are just file names; data directory needs to be entered in `config.yaml`

`sra_result.csv` (md5: 9b99cbddaf8bda4c220cba016c620db6) - metadata obtained from NCBI SRA site. From [list of 58 SRA experiments](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=542925) choose: Send to - File - Summary

`SraAccList.txt` (md5: dad41b5873fe3e85bf1e658f804ab54d) - list of just the SRA accession numbers, potentially useful for downloading files. Obtained from NCBI SRA site. From [list of 58 SRA experiments](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=542925) choose: Send to - File - Accession list

`SraRunInfo.csv` (md5: 76612079dc3afe475a00fde2df04fb68) - metadata obtained from NCBI SRA site. From [list of 58 SRA experiments](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=542925) choose: Send to - File - RunInfo
