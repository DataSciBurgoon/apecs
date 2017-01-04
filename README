# Autoencoder Predicting Estrogenic Chemical Substances (APECS)

# Author and Other Misc Info
Lyle D. Burgoon, Ph.D. US Army Engineer Research and Development Center. Research Triangle Park, NC.

# License
Public domain -- this is a work of the United States Government.

# Pre-Reqs
1. You must have R installed (at least version 3.2.4)
2. Java (64-bit) version 1.6 or later

# Usage
The easiest way to make predictions is to name your tab-delimited input file (see below for more information about the input file) data.txt, and put it in the in\_vitro\_script directory (to use the model trained on in vitro calls) or the in\_vivo\_script directory (to use the model trained on the in vivo calls). Then you can just double click apecs\_in\_vitro.bat or apecs\_in\_vivo.bat (on Windows) or run apecs\_in\_vitro.sh or apecs\_in\_vivo.sh on the command line (Linux/Mac).

More advanced users can run RScript on the command line to run either the apecs\_in\_vitro.R or the apecs\_in\_vivo.R scripts. Or, they could use R directly or RStudio to run the R scripts.

# Input File
The demo_data.txt file shows the input we're expecting. This was taken directly from the ToxCast MySQL database. The R scripts are expecting a tab-delimited input file. The columns are:

1. chnm: chemical name
2. casn: CAS registration number
3. gene_name: the gene name from ToxCast
4. assay_name: the assay name from ToxCast
5. assay_format: the format for the assay
6. assay_source: the source for the assay
7. organism: the organism that the assay was derived from
8. logc: the log of the chemical's concentration for that row
9. resp: the response value
10. idx: a remnant column that I added for a different purpose -- you can ignore this column

Note that chemicals are repeated in the rows. This is because each row represents the response value for a given chemical at a given concentration in a given assay. 

# Do I Have To Use ToxCast Data?
Nope. Frankly, the software doesn't care where your data comes from. All that matters is that it's formatted the same way as the demo file. You'll notice that several of these columns aren't actually used by the software. The file format you're seeing is the default format that we use for all of our ToxCast database queries.



