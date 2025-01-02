import pandas as pd

# Load your data
oneQ = pd.read_csv("path/to/your/file/deseq2_Ksg1qvsNDCresults.csv")
threeQ = pd.read_csv("path/to/your/file//deseq2_Ksg3qvsNDCresults.csv")
frank = pd.read_csv("path/to/your/file/5seqBeds/KsgData_fRanked.csv")
rnaseq = pd.read_csv("path/to/your/RNAseqResultsFile/T4_Ksg_RNAseq_1Q.csv")

# Initialize new columns
frank["fiveseqDEseq2"] = ""
frank["fiveseqpadj"] = ""
frank["RNAseqDEseq2"] = ""
frank["RNAseqpadj"] = ""
frank["location"] = ""  # Initialize location column
frank.insert(frank.columns.get_loc("TSSUTR"), "paste", frank["TSSUTR"].astype(str) + frank["Locus"].astype(str))

for i in range(len(frank)):
    # Check against oneQ data frame
    if frank.loc[i, "Genome"] in ["Ksg1q_t60", "NDCt60"]:
        match_oneQ = oneQ[oneQ.iloc[:, 0] == frank.loc[i, "SuperPos"]]
        if not match_oneQ.empty:
            frank.loc[i, "fiveseqDEseq2"] = match_oneQ["log2FoldChange"].values[0]
            frank.loc[i, "fiveseqpadj"] = match_oneQ["padj"].values[0]
    
    # Check against threeQ data frame
    if frank.loc[i, "Genome"] == "Ksg3q_t60":
        match_threeQ = threeQ[threeQ.iloc[:, 0] == frank.loc[i, "SuperPos"]]
        if not match_threeQ.empty:
            frank.loc[i, "fiveseqDEseq2"] = match_threeQ["log2FoldChange"].values[0]
            frank.loc[i, "fiveseqpadj"] = match_threeQ["padj"].values[0]
    
    # Check against RNAseq DataFrame
    match_RNAseq = rnaseq[rnaseq["Gene"] == frank.loc[i, "paste"]]
    if not match_RNAseq.empty:
        frank.loc[i, "RNAseqDEseq2"] = match_RNAseq["Value"].values[0]  # Ensure correct column name
        frank.loc[i, "RNAseqpadj"] = match_RNAseq["padj"].values[0]
    
    # Determine location based on TSSUTR and Locus
    if pd.notna(frank.loc[i, "TSSUTR"]) and frank.loc[i, "TSSUTR"]:
        frank.loc[i, "location"] = "5'UTR"
    elif pd.notna(frank.loc[i, "Locus"]) and frank.loc[i, "Locus"]:
        frank.loc[i, "location"] = "Genic"

# Delete index 1 (paste column)
frank.drop(frank.columns[1], axis=1, inplace=True)

# Write the resulting DataFrame to a CSV file
frank.to_csv("/Users/merinakzo/Desktop/Ksg_fRank_DEseq2_python.csv", index=False)
