To run: open "Analyze_barcodes.ipynb", download the appropriate packages (Bio, pandas, os, numpy, collections). Then input the directory containing your files where it says "Directory_here" and run.
The directory should contain the following files and folders:
1. A "Fastq_folders" folder with fastq files (eg fname.fastq.gz)
2. "Barcode_info.csv" with a "Barcode_ID" column (with the names of each barcode) and a "Barcode_seq" column (with the actual sequence of each barcode)
3. "Fastq_file_IDs.csv" with a "Filename" column (e.g. fname.fastq, leave off the ".gz") and a "Sample_label" column (whatever name you have for that sample)
4. "Fastq_file_IDs_UMI.csv" that's exactly the same as Fastq_file_IDs.csv
5. "Flanks.csv" with cell A1 = "Left_flank", cell B1 = "Right_flank", and A2 and B2 the left and right sequence flanks for the barcode respectively. e.g. if the barcode (call it NNNN) is located in the sequence context "ACGTTTCNNNNCGCACCG", then left flank is ACGTTTC and right flank is CGCACCG. Make the flank long enough to uniquely specify a region of the amplicon.
6. "Read_tags.csv" that has column "Read_tag" and column "Tag_descriptor". Read_tag is a tag that I insert into the primers between the sequencing adapter (for Azenta Amplicon EZ sequencing) and the annealing portion of the primer that allows for multiplexing, i.e. the ACGT read tag primer would look like "Adapter-ACGT-Primer". Tag_descriptor is whatever name you want to give.
7. "UMI_flanks.csv": the flanks of the UMI sequence, structured the same as Flanks.csv but for where the UMI is inserted.
