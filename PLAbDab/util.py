import pandas as pd


def add_url_to_paired_data(df):
    links = []

    for entry in df.iloc:

        if "patent US" in entry.heavy_definition:
            patent_id = entry.heavy_definition.split()[-1]
            link = "https://patents.google.com/patent/US{}/en".format(patent_id)
        elif entry.pairing == "Xtal structure":
            link = "https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/structureviewer/?pdb={}".format(entry.ID[:4].lower())
        elif entry.pairing == "TheraSAbDab":
            link = "https://opig.stats.ox.ac.uk/webapps/newsabdab/therasabdab/search/?therapeutic={}".format(entry.ID.replace("_2",""))
        elif entry.pairing == "Same entry":
            link = "https://www.ncbi.nlm.nih.gov/protein/{}".format(entry.ID)
        else:
            link = "https://www.ncbi.nlm.nih.gov/protein/{}".format(entry.heavy_ID)

        links.append(link)

    df = df.copy()
    df["url"] = links
    return df