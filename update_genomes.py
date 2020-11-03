#! /usr/bin/python

import os, urllib.request, requests
from datetime import date
from multiprocessing.pool import ThreadPool
import sys

print(sys.getdefaultencoding())

""" Globals variables"""
# Folder to store the complete GENOME_REPORT containing all the genomes available at NCBI\
# and the resumed summary with the genomes included in the Genome collection
tmp = 'Genomes/tmp'
genome_dir = 'Genomes/Bacteria'
url_prok = 'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt'
taxons = 'taxons.txt'


# << ACCESSORY FUNTIONS >>
# Replace special characters at species names
def replace(str):

	replace = ["-", ";", " ", "|", ".", "\\", ":", ",", "'", 
               "&", "~", "@", "+", "^", "°", "$", "=", "*", 
               "!", "µ", "(", ")", "[", "]", "__", "___", "\#"]
	for x in replace:
		if x in str:
			str = str.replace(x, "_")
	return str


# Checks if the URL requested exists
def url_exists(path):
	r = requests.head(path)
	return r.status_code == requests.codes.ok


# Download URLs if files don't exist.
def fetch_url(entry):
	path, url = entry
    if os.path.exists(path):
        print("\tFile",path,"already exists.")
        pass
    else:
		r = requests.get(url, stream=True)
		if r.status_code == 200:
			with open(path, 'wb') as f:
				for file in r:
					f.write(file)
		
	return path


# << MAIN FUNCTIONS >>
# Create output folders: Genomes and Temp
def create_folders(dir_list):

    for folder in dir_list:
        try:
            makedirs(folder)
        except:
            print("The output directory %s already exists." % folder)
            pass


# Download the Prokaryotes summary with all genomes available at NCBI repository
def download_summary(url, tmp_dir, change_name = None):

	summary_name = str(date.today())+ "_" + url.split("/")[-1]
	summary_path = "/".join([tmp_dir, summary_name])
    
	if os.path.exists(summary_path):
		print("This file already exists.", os.path.realpath(summary_path))
	else:
		urllib.request.urlretrieve(url, summary_path)

	return summary_path


# Read file containing the names of Genus/Species to be downloaded
def get_taxons(table, report, tmp_dir):
	""" Headers of taxa table [0]Superkingdom [1]Kingdom [2]Phylum [3]Class [4]Order [5]Family [6]Genus [7]Species """
	list_taxa = []
	with open(table, "r") as handle:
		next(handle)
		for line in handle:

			if line.startswith("#"):
				pass #continue
			else:
				data = line.strip().split("\t")
				genus = data[6]
				species = data[7]
				# Species = '*' means that all species from genus will be downloaded
				if species == "*":
					list_taxa.append(genus)
				# Downloaded genomes at species level
				else:
					name = " ".join([genus, species])
					list_taxa.append(name)
	print(list_taxa, "\n")

	return list_taxa 


def select_from_list(tax_list, report, tmp_dir):

	report_path = os.path.join(tmp_dir, str(date.today())+ "_" + "downloaded.txt")
	dload_report = open(report_path, "w")
	subject = open((report), "r", encoding="UTF8")
	headers = subject.readline().strip() + "\tUrlGCA\tNameGCA\tUrlGCF\tNameGCF\tUrlWGS\tNameWGS\n"
	dload_report.write(headers)

	""" Headers/Fields in prokaryotes.txt
	[0]Organism/Name [1]TaxID [2]BioProjectAccession [3]BioProjectID [4]Group [5]SubGroup [6]Size(Mb)
	[7]GC% [8]Replicons [9]WGS [10]Scaffolds [11]Genes [12]Proteins [13]ReleaseDate [14]ModifyDate
	[15]Status [16]Center [17]BioSampleAccession [18]Assembly Accession [19]Reference [20]FTP Path
	[21]PubmedID [22]Strain + {Columns added after 'get_taxons'} [23]UrlGCA [24]NameGCA [25]UrlGCF [26]NameGCF [27]UrlWGS_GCA [28]NameWGS_GCA
	"""

	for line in subject:
		for t in tax_list:
			if t in str(line.strip()):

				fields = line.strip().split("\t")

				# Use "https" to download file instead of "ftp" bercause the HCP doesn't accepts it.
				url = fields[20].replace("ftp:", "https:")
				file_GCA = url + "/" + url.split("/")[-1] + "_genomic.gbff.gz"
				file_WGS_GCA = url + "/" + url.split("/")[-1] + "_wgsmaster.gbff.gz"

				# For GCF files
				url_GCF = url.replace("GCA", "GCF")
				file_GCF = url_GCF + "/" + url_GCF.split("/")[-1] + "_genomic.gbff.gz"
				
                name = fields[0].replace("uncultured ", "")
				genus = replace(name.split(" ")[0]).strip()
				label = genus[0]
				sp = name.split(" ")[1]
				wgs = fields[9] # accession number of WGS

				if fields[22] != "-":
					strain = fields[22]
				elif strain == "-":
					if wgs != "-":
						strain = wgs

				strain = replace(strain.strip())

				if "subsp" in name:  # Include 'subsp.' in the name
					# Output file name
					subsp = name.split(" ")[3]

					file_name_GCA = "_".join([genus, sp, "subsp.", subsp, strain, "GCA.gbff.gz"])
					file_path_GCA = "Genomes/Bacteria/" + genus + "/" + str("_".join([label,sp,subsp]) + "/" + strain + "/" + file_name_GCA)

					file_name_GCF = "_".join([genus, sp, "subsp.", subsp, strain, "GCF.gbff.gz"])
					file_path_GCF = "Genomes/Bacteria/" + genus + "/" + str("_".join([label,sp,subsp]) + "/" + strain + "/" + file_name_GCF)

					file_name_WGS = "_".join([genus, sp, "subsp.", subsp, strain, "WGS.gbff.gz"])
					file_path_WGS = "Genomes/Bacteria/" + genus + "/" + str("_".join([label,sp,subsp]) + "/" + strain + "/" + file_name_WGS)

					dload_report.write(line.strip() + "\t" + file_GCA + "\t" + file_path_GCA + "\t" + file_GCF + "\t" + file_path_GCF + "\t" + file_WGS_GCA + "\t" + file_path_WGS + "\n")
				else:
					file_name_GCA = "_".join([genus, sp, strain, "GCA.gbff.gz"])
					file_path_GCA = "Genomes/Bacteria/" + genus + "/" + str("_".join([label,sp]) + "/" + strain + "/" + file_name_GCA)

					file_name_GCF = "_".join([genus, sp, strain, "GCF.gbff.gz"])
					file_path_GCF = "Genomes/Bacteria/" + genus + "/" + str("_".join([label,sp]) + "/" + strain + "/" + file_name_GCF)

					file_name_WGS = "_".join([genus, sp, strain, "WGS.gbff.gz"])
					file_path_WGS = "Genomes/Bacteria/" + genus + "/" + str("_".join([label,sp]) + "/" + strain + "/" + file_name_WGS)

					dload_report.write(line.strip() + "\t" + file_GCA + "\t" + file_path_GCA + "\t" + file_GCF + "\t" + file_path_GCF + "\t" + file_WGS_GCA + "\t" + file_path_WGS + "\n")

	return report_path


#TODO compare_summaries()
""" Here the scripts will compare the new summary and the latest downloaded to detect the new genomes included """



# tbl = downloaded.txt - Table with filtered genomes based on the taxa listed in 'taxons.txt'.
def get_genomes(tbl):
	with open(tbl, "r") as handle:

		for line in handle:

			if line.startswith("#"):
				pass
			else:
				fields = line.strip().split("\t")
				out_dir = "/".join(fields[24].split("/")[:-1])

				try:
        				os.makedirs(out_dir)
				except:
					print("#", out_dir, "already exists")
                
                #store the URLs of genomes that will be download
				urls = []

				url_GCA, file_path_GCA = fields[23], fields[24]
				urls.append((file_path_GCA, url_GCA))

				url_GCF, file_path_GCF = fields[25], fields[26]
				urls.append((file_path_GCF, url_GCF))

				url_WGS_GCA, file_path_WGS_GCA = fields[27], fields[28]
				urls.append((file_path_WGS_GCA, url_WGS_GCA))

				ThreadPool(5).imap_unordered(fetch_url, urls)


# Main steps
create_folders([tmp, genome_dir])
summary = download_summary(url_prok, tmp)
dload = get_taxons(taxons, summary, tmp)
selected = select_from_list(dload, summary, tmp)
get_genomes(selected)
