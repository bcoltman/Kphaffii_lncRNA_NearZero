"""
Copyright [2009-present] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
Usage: python search.py [file] [database]
Examles:
python search.py file.fasta  # (search in the RNAcentral database)
python search.py file.fasta mirbase
"""

import json
import math
import requests  # pip install requests
import sys
import time
from Bio import SeqIO  # pip install biopython
from pathlib import Path
import pandas as pd

# sequence search server
server = "https://search.rnacentral.org/"



def get_sequence_search_result(description, job_id, e_value, percent_identity):
	"""
	Function to check the status of the job and also to get the results of the sequence search
	:param description: description of the FASTA record
	:param sequence: sequence of the FASTA record
	:param job_id: id of the job
	:return: dict containing the results
	"""
	get_status = requests.get(server + "api/job-status/" + job_id)
	job_status = json.loads(get_status.text)["status"]
	if job_status == "success" or job_status == "partial_success":
		# get results
		# there is a limit on the number of results that can be requested
		start = 0
		get_result = requests.get(server + "api/facets-search/" + job_id + "?start=" + str(start) + "&size=100")
		if get_result.status_code == 502:
			print("Bad Gateway")
			print("There was an error in the following record: {}".format(description))
			return {"job_id": job_id, "status": job_status, "description": description}
		else:
			get_result = json.loads(get_result.text)
			results = [get_result["entries"]]
			# check the number of similar sequences and make new requests if necessary
			hit_count = get_result["hitCount"]
			iter_number = int(math.ceil(hit_count / 100.0))
			for num in range(iter_number - 1):
				start += 100
				new_request = requests.get(server + "api/facets-search/" + job_id + "?start=" + str(start) + "&size=100")
				if new_request and new_request.text:
					new_request_result = json.loads(new_request.text)
					results.append(new_request_result["entries"])
			results = [item for sublist in results for item in sublist] 			
			#results = [item for sublist in results for item in sublist if (item["e_value"] < e_value)&(item["identity"] > percent_identity)] #print(results[0][0]["fields"]["tax_string"])
		job_result = False
		if results:
			job_result = {"job_id": job_id,"status": job_status,"description": description,"results": results}
		return job_result
	elif job_status == "error" or job_status == "timeout":
		# return the metadata of the search
		# you can also send us the job_id so we can investigate further
		print("There was an error in the following record: {}".format(description))
		return {"job_id": job_id, "status": job_status, "description": description}
	elif job_status == "started" or job_status == "pending":
		# try again in 10 seconds
		time.sleep(10)
		return get_sequence_search_result(description, job_id, e_value, percent_identity)


def get_rfam_result(job_id, e_value):
	"""
	Function to check the status of the Rfam job and also to get the results
	:param job_id: id of the job
	:return: rfam results
	"""
	infernal_status = requests.get(server + "api/infernal-status/" + job_id)
	infernal_status = json.loads(infernal_status.text)["status"]
	results = []
	if infernal_status == "success":
		# get results
		infernal_results = requests.get(server + "api/infernal-result/" + job_id)
		infernal_results = json.loads(infernal_results.text)
		if infernal_results:
			for item in infernal_results:
				
				results.append({
				    "family": item["description"],
				    "accession": item["accession_rfam"],
				    "start": item["seq_from"],
				    "end": item["seq_to"],
				    "bit_score": item["score"],
				    "e_value": item["e_value"],
				    "strand": item["strand"],
				    "alignment": item["alignment"]
				})
				#if item["e_value"] < e_value:
				#	results.append({
				#	    "family": item["description"],
				#	    "accession": item["accession_rfam"],
				#	    "start": item["seq_from"],
				#	    "end": item["seq_to"],
				#	    "bit_score": item["score"],
				#	    "e_value": item["e_value"],
				#	    "strand": item["strand"],
				#	    "alignment": item["alignment"]
				#	})
			return results if results else False
	elif infernal_status == "error" or infernal_status == "timeout":
		return "There was an error in the Rfam search"
	elif infernal_status == "started" or infernal_status == "pending":
		# try again in 10 seconds
		time.sleep(10)
		return get_rfam_result(job_id, e_value)



def main():
	filename = None
	database = []
	e_value = 1e-5
	percent_identity = 90


	if len(sys.argv) != 5:
		print("You must specify filename, directory, e_value and percent_identity threshold")
		exit()

	elif len(sys.argv) == 5:
		filename = sys.argv[1]
		directory = sys.argv[2]
		e_value = float(sys.argv[3])
		percent_identity = float(sys.argv[4])

	else:
		print("Usage: python search.py file.fasta")
		exit()


	sequence_df = pd.DataFrame(columns=["Transcript","Gene","ID","Source","RNA_type","Organism","Query Coverage", "Target Coverage", "E-value","Score","Identity","Taxonomy"])
	rfam_df = pd.DataFrame(columns=["Transcript","Gene", "Family","accession","start","end","bit_score","e_value","strand","alignment"])

#	submission_df = pd.DataFrame(columns=["JobID", "Description", "Sequence"])
#
#	with open(filename, mode='r') as handle:
#	# use Biopython's parse function to process individual FASTA records
#		for record in SeqIO.parse(handle, 'fasta'):
#			# extract individual parts of the FASTA record
#			description = record.description
#			transcript, gene = description.split()
#			sequence = record.seq
#			# submit a job
#			data = {"databases": database, "query": str(sequence)}
#			post_job = requests.post(server + "api/submit-job", json=data)
#			# get job_id
#			job_id = None
#			if post_job.status_code == 201:
#				job_id = json.loads(post_job.text)["job_id"]				
#			else:
#				time.sleep(10)
#				post_job = requests.post(server + "api/submit-job", json=data)
#				if post_job.status_code == 201:
#					job_id = json.loads(post_job.text)["job_id"]
#				else:
#					print("Failed to submit job. Record that failed:\n {} \n {}".format(description, sequence))
#			submission_df = pd.concat([submission_df, pd.DataFrame.from_dict([{"JobID":job_id, "Description":description, "Sequence":str(sequence)}])])
#
#
#	submission_df.to_csv(f"{directory}/RNACentralJobs.tsv", sep = "\t", index=False)
	submission_df = pd.read_csv(f"{directory}/RNACentralJobs.tsv", sep = "\t")

	for index, row in submission_df.iterrows():

		# Extract the data from the first column (index 0)
		job_id, description, sequence = row
		transcript, gene = description.split()
		if job_id:
					# get sequence search results
			sequence_search = get_sequence_search_result(description, job_id, e_value, percent_identity)
		if sequence_search:		
			sequence_df = pd.concat([sequence_df, pd.DataFrame.from_dict([{"Transcript":transcript, 
							"Gene":gene,
							"ID":hit["id"], 
							"Source":hit["source"], 
							"RNA_type":hit["fields"]["rna_type"], 
							"Organism":hit["fields"]["standard_name"], 
"Query Coverage":hit["query_coverage"], "Target Coverage":hit["target_coverage"],
							"E-value":hit["e_value"], 
							"Score":hit["score"], 
							"Identity":hit["identity"], 
							"Taxonomy":hit["fields"]["tax_string"]} for hit in sequence_search["results"]])]) #"Alignment":hit["alignment"]
		if job_id:# get Rfam results
			rfam = get_rfam_result(job_id, e_value)
		if rfam:
			rfam_df = pd.concat([rfam_df, pd.DataFrame.from_dict([{"Transcript":transcript, 
					"Gene":gene,
					"Family":hit["family"], 
					"accession":hit["accession"], 
					"start": hit["start"],
					"end": hit["end"],
					"bit_score": hit["bit_score"],
					"e_value": hit["e_value"],
					"strand": hit["strand"],
					"alignment": hit["alignment"]} for hit in rfam])])
		print(f"Done - {description}")


	sequence_df.to_csv(f"{directory}/RNACentralAll.tsv", sep = "\t", index=False)
	rfam_df.to_csv(f"{directory}/Rfam_FamiliesAll.tsv", sep = "\t", index=False)
	
	sequence_df.loc[(sequence_df["E-value"] < e_value)&(sequence_df["Identity"] > percent_identity),:].to_csv(f"{directory}/RNACentral.tsv", sep = "\t", index=False)
	rfam_df.loc[(rfam_df["e_value"] < e_value),:].to_csv(f"{directory}/Rfam_Families.tsv", sep = "\t", index=False)


if __name__ == "__main__":
	main()
